#ifndef bdpt_h__
#define bdpt_h__

#include <vector>
#include "camera.h"
#include "consoledebug.h"
#include "film.h"
#include "geometry.h"
#include "intersection.h"
#include "lightpath.h"
#include "material.h"
#include "numeric.h"
#include "sample.h"
#include "scene.h"
#include "triangle.h"

// reference BDPT implementation based on Hachisuka's MMLT demo
// this implementation is straightforward but less efficient
class ReferenceBDPT
{
    const Camera&            camera;
    FrameBuffer&             film;
    const Scene&             scene;
    std::vector<FrameBuffer> pyramid;

private:
    // path sampling (recursively construct a path, but the first vertex is
    // always constructed outside this function)
    void TracePath(Path&            path,
                   const Ray&       ray,
                   const int        RayLevel,
                   const int        MaxRayLevel,
                   PathType         path_type,
                   StandardSampler& sampler)
    {
        if (RayLevel >= MaxRayLevel)
        {
            return;
        }

        HitInfo hit = scene.intersect(ray);
        if (!hit)
        {
            return;  // to background
        }

        vec3 position = ray.at(hit.distance());
        vec3 object_normal = hit.shadingNormal();
        vec3 orienting_normal =
            dot(object_normal, ray.dir) < 0 ? object_normal : -object_normal;

        const Vert& last_vertex = path.end();
        float       total_pdfA = last_vertex.getTotalPdfA();
        vec3        throughput = last_vertex.getThroughput();

        if (RayLevel == 1 && path_type == PathType::FromEye)
        {
            total_pdfA *=
                camera.measurementFunction(position, orienting_normal);
        }
        else
        {
            const vec3& previous_position =
                last_vertex.p;  // can also use ray.o
            float previous_pdfW = last_vertex.pWo;
            float pdfA = previous_pdfW *
                         pdfWtoA(previous_position, position, orienting_normal);
            total_pdfA *= pdfA;

            float G = geometryTerm(previous_position,
                                   last_vertex.getOrientingNormal(),
                                   position,
                                   orienting_normal);
            throughput *= last_vertex.weight * G;
        }

        const TriangleObject* obj = hit.triangleObject();
        const Material&       mat = todo_getMaterial(obj->materialID());

        const float rnd0 = sampler.next();
        const float rnd1 = sampler.next();

        vec3  total_weight;
        vec3  bsdf_weight;
        float pdfW;
        float cos_wo;

        Ray nr;
        nr.orig = position;
        // sample a reflected direction according to the scattering material
        nr.dir = mat.sample(-ray.dir,
                            object_normal,
                            rnd0,
                            rnd1,
                            total_weight,
                            &pdfW,
                            &cos_wo,
                            &bsdf_weight);
        if (mat.isDelta())
        {
            bsdf_weight /= cos_wo;
        }

        // new orientation in case refraction
        orienting_normal =
            dot(object_normal, nr.dir) > 0 ? object_normal : -object_normal;

        // set path data
        Vert vert(position, object_normal, obj);
        {
            vert.pWo = pdfW;
            vert.pWi = sqrt(-1);  // TODO : use static_NAN
            vert.weight = bsdf_weight;
            vert.cos_wo = cos_wo;
            vert.debug_wi = -ray.dir;
            vert.debug_wo = nr.dir;
            vert.debug_type = path_type;
            vert.setDeltaFlag(mat.isDelta());
            vert.setTotalPdfA(total_pdfA);
            vert.setOrientingNormal(orienting_normal);
            vert.setThroughput(throughput);
            vert.debug_vertex_type =
                mat.isLight() ? Emitter : (mat.isDelta() ? Specular : Diffuse);
        }

        path.insert(vert);

        // recursive ray tracing for scattering surfaces
        // TODO : transform to iterative
        if (!mat.isLight())
        {
            TracePath(path, nr, RayLevel + 1, MaxRayLevel, path_type, sampler);
        }
    }

    Path GenerateEyePath(const int MaxEyeEvents, StandardSampler& sampler)
    {
        Path Result;

        if (MaxEyeEvents == 0)
        {
            return Result;
        }

        Ray r = camera.makeRay(sampler.next(), sampler.next());

        Vert v(r.orig, camera.forwardDirection(), nullptr);
        v.debug_type = FromEye;
        v.debug_wo = r.dir;
        v.debug_vertex_type = Lens;
        v.setDeltaFlag(false);
        v.setTotalPdfA(1.0);
        v.setOrientingNormal(camera.forwardDirection());
        v.setThroughput(vec3(1.0));
        v.weight = vec3(1);
        v.pWo = 1.0;
        Result.insert(v);  // first vertex is on the camera aperture (if pinhole
                           // then it's the camera position)

        TracePath(Result, r, 1, MaxEyeEvents, FromEye, sampler);
        return Result;
    }

    Path GenerateLightPath(const int MaxLightEvents, StandardSampler& sampler)
    {
        Path Result;

        if (MaxLightEvents == 0)
        {
            return Result;
        }

        const TriangleObject* light_obj = nullptr;
        vec3                  position, normal, direction;
        float                 pdfA, pdfW;
        scene.sampleLightAreaAndDirection(
            position, normal, pdfA, direction, pdfW, light_obj, sampler);
        assert(nullptr != light_obj);

        Vert v(position, normal, light_obj);
        v.debug_type = FromLight;
        v.debug_wo = direction;
        v.debug_vertex_type = Emitter;
        v.setDeltaFlag(false);
        v.setTotalPdfA(pdfA);
        v.setOrientingNormal(normal);
        v.setThroughput(todo_getMaterial(light_obj->materialID())
                            .evaluate(vec3(), normal, direction));
        v.weight = vec3(1);
        v.pWo = pdfW;
        Result.insert(v);  // first vertex is on the light's surface

        Ray r(position, direction);
        TracePath(Result, r, 1, MaxLightEvents, FromLight, sampler);
        return Result;
    }

    // check if the path can be connected or not (visibility term)
    // note that there would be false positives if eps is not set,
    // and this would lead to fireflies for bsdf may return a zero valued cosine
    // term
    bool isConnectable(const Path& eye, const Path& light, float& px, float& py)
    {
        vec3        Direction;  // primary ray direction for rasterization
        const Vert& ev =
            eye.end();  // the last vertex of eye path, to be connected
        const Vert& lv =
            light.end();  // the last vertex of light path, to be connected

        // cannot be connected
        if (ev.isDelta() || lv.isDelta())
        {
            return false;
        }

        bool be_visible;

        // no direct hit to the film (pinhole)
        if ((eye.size() == 0) && (light.size() >= 2))
        {
            be_visible = false;
        }

        // direct hit to the light source
        // ( the only ways to directly render the light are
        // Xeye.n == 2 && Xlight.n == 0, and Xeye.n == 1 && Xlight.n == 1 )
        else if ((eye.size() >= 2) && (light.size() == 0))
        {
            be_visible = todo_getMaterial(ev.obj->materialID()).isLight();

            // a temporary fix for the fireflies due to s=0 paths that have
            // second last vertices too close to the light
            // TODO : see if s=0 paths are difficult for BDPT
            //             be_visible =
            //             todo_getMaterial(ev.obj->materialID()).isLight() &&
            //             ev.canConnectTo(eye[eye.size() - 2].p);

            Direction = normalize(eye[1].p - eye[0].p);
        }

        // direct connection to camera
        else if ((eye.size() == 1) && (light.size() >= 1))
        {
            Ray r(eye[0].p, normalize(lv.p - eye[0].p));
            be_visible = !scene.occluded(r, length(eye[0].p - lv.p));
            Direction = r.dir;
        }

        // shadow ray connection
        else
        {
            vec3 d = normalize(lv.p - ev.p);
            Ray  r(ev.p, d);
            be_visible = ev.canConnectTo(lv.p) && lv.canConnectTo(ev.p) &&
                         !scene.occluded(r, length(ev.p - lv.p));
            Direction =
                normalize(eye[1].p - eye[0].p);  // primary ray direction
        }

        // get the pixel location
        return be_visible &&
               camera.rasterizePrimaryRay(
                   Direction, px, py);  // the path must be rasterizable
    }

    // combine two traced paths into all possible subpaths and return pixel
    // contributions
    /// BPT connections (loop over the path pyramid)
    /// - limit the connection to a specific technique if s and t are provided
    PathContribution CombinePaths(const Path& EyePath,
                                  const Path& LightPath,
                                  const int   SpecifiedNumEyeVertices = -1,
                                  const int   SpecifiedNumLightVertices = -1)
    {
        PathContribution Result;
        Result.n = 0;
        Result.sc = 0.0;
        const bool Specified = (SpecifiedNumEyeVertices != -1) &&
                               (SpecifiedNumLightVertices != -1);

        // MaxEvents = the maximum number of vertices
        for (int PathLength = MinPathLength; PathLength <= MaxPathLength;
             PathLength++)
        {
            std::vector<double> debug_weights;

            for (int NumEyeVertices = 0; NumEyeVertices <= PathLength + 1;
                 NumEyeVertices++)  // note that num verts = path length + 1
            {
                const int NumLightVertices = (PathLength + 1) - NumEyeVertices;

                if (NumEyeVertices == 0)
                    continue;  // no direct hit to the camera (pinhole camera,
                               // aperture size = 0)
                if (NumEyeVertices > EyePath.size())
                    continue;  // impossible case
                if (NumLightVertices > LightPath.size())
                    continue;  // impossible case

                // take only the specified technique if provided
                if (Specified &&
                    ((SpecifiedNumEyeVertices != NumEyeVertices) ||
                     (SpecifiedNumLightVertices != NumLightVertices)))
                {
                    continue;
                }

                // extract subpaths with specific length
                Path Eyesubpath = EyePath.subpath(NumEyeVertices);
                Path Lightsubpath = LightPath.subpath(NumLightVertices);

                // check the path visibility
                float px = -1.0,
                      py = -1.0;  // pixel contribution position in pixel units
                if (!isConnectable(Eyesubpath, Lightsubpath, px, py))
                {
                    continue;
                }

                // construct a full path (Eyesubpath + Lightsubpath ->
                // SampledPath)
                CombinedPath SampledPath(Eyesubpath, Lightsubpath);

                //////////////////////////////////////////////////////////////////////////
                // main evaluation
                // evaluate the path
                vec3d f = vec3d(pathThroughputWeighted(SampledPath));

                // this could be tricky
                //                 vec3 ff =
                //                 SampledPath.eyeEnd().getThroughput(); if
                //                 (SampledPath.lightSize() == 0)
                //                 {
                //                     ff *=
                //                     SampledPath.lightEnd().getThroughput();
                //                 }
                //                 else
                //                 {
                //                     ff *=
                //                     todo_getMaterial(SampledPath.eyeEnd()->getMaterialID()).getBSDF()->evaluate(vec3(),
                //                     normal, vec3())
                //                 }
                //                 vec3d f = vec3d(ff);

                // uniform weighting
                // this also seems to work
                // int num_non_specs = SampledPath.numNonSpecularEdges();
                // double path_weight = (1.0 / (num_non_specs + 1));

                // path mis
                double path_weight = calculateMisWeight(
                    SampledPath);  // double precision is important
                //                 assert(path_weight == path_weight);
                if (path_weight != path_weight)
                {
                    path_weight = 0;
                }

                vec3 c = (f * path_weight)
                             .toFloat();  // Monte Carlo integration formula
                //////////////////////////////////////////////////////////////////////////

                if (isZero(c))
                {
                    continue;
                }

                // store the pixel contribution (pushback style)
                Result.c[Result.n] =
                    Contribution(px, py, c, NumLightVertices, NumEyeVertices);
                Result.n++;

                // scalar contribution function
                Result.sc = f_max(maxComponent(c), Result.sc);

                // return immediately if the technique is specified
                if (Specified && (SpecifiedNumEyeVertices == NumEyeVertices) &&
                    (SpecifiedNumLightVertices == NumLightVertices))
                {
                    return Result;
                }
            }
        }
        return Result;
    }

    // accumulate the pixel contributions of a path onto the film (note that
    // normalization for sample numbers is done when writing to image)
    void AccumulatePathContribution(const PathContribution pc)
    {
        if (pc.sc == 0)
        {
            return;
        }

        for (int i = 0; i < pc.n; i++)
        {
            const int ix = int(pc.c[i].x),
                      iy = int(pc.c[i].y);  // pixel index, no advanced
                                            // reconstruction filter is used
            const vec3 c = pc.c[i].c;       // scaled pixel contribution

            film.accumulatePixel(ix, iy, c);
        }

        for (int i = 0; i < pc.n; i++)
        {
            const int ix = int(pc.c[i].x),
                      iy = int(pc.c[i].y);  // pixel index, no advanced
                                            // reconstruction filter is used
            const vec3 c = pc.c[i].c;       // scaled pixel contribution

            int s = pc.c[i].s;
            int t = pc.c[i].t;

            pyramid[s * MaxEvents + t].accumulatePixel(ix, iy, c);
        }
    }

    double calculatePdfA(const std::vector<const Vert*>& vs,
                         const int                       from_idx,
                         const int                       next_idx)
    {
        const Vert& from_vertex = *vs[from_idx];
        const Vert& next_vertex = *vs[next_idx];
        const int   prev_from_idx = (from_idx - next_idx) + from_idx;
        Vert const* prev_from_vertex = nullptr;
        if (0 <= prev_from_idx && prev_from_idx < vs.size())
        {
            prev_from_vertex = vs[prev_from_idx];
        }

        const vec3 to = next_vertex.p - from_vertex.p;
        const vec3 normalized_to = normalize(to);
        double     pdf = 0;

        switch (from_vertex.debug_vertex_type)
        {
            default:
            case Lens:
                return camera.measurementFunction(
                    next_vertex);  // pdfA of sampling the next vertex
                break;
            case Emitter:
            case Diffuse:
                pdf = todo_getMaterial(from_vertex.obj->materialID())
                          .pdfW(vec3(0, 0, 0), from_vertex.n, normalized_to);
                break;
            case Specular:
                //             pdf = 1;
                pdf = from_vertex.pWo;
                break;
        }

        // pdfW to pdfA
        const vec3 next_new_orienting_normal =
            dot(to, next_vertex.n) < 0 ? next_vertex.n : -next_vertex.n;
        return pdf *
               (dot(-normalized_to, next_new_orienting_normal) / dot(to, to));
    }

    double calculateMisWeight(const CombinedPath& path)
    {
        int t = path.eyeSize();
        int s = path.lightSize();

        std::vector<double>      pi1_pi(t + s);  // p_i+1 over p_i
        std::vector<const Vert*> vs(t + s);
        double                   PA_y0 = 1.0 / scene.lightArea();
        double                   PA_x0 = 1.0;

        std::vector<std::pair<double, double> > debug_p(s + t);

        const int k = t + s - 1;
        for (int i = 0; i < path.size(); i++)
        {
            vs[i] = &path[k - i];
        }

        pi1_pi[0] = PA_y0 / (calculatePdfA(vs, 1, 0));
        debug_p[0] = std::make_pair(PA_y0, calculatePdfA(vs, 1, 0));
        for (int i = 1; i < k; i++)
        {
            double a = calculatePdfA(vs, i - 1, i);
            double b = calculatePdfA(vs, i + 1, i);
            pi1_pi[i] = a / b;
            debug_p[i] = std::make_pair(a, b);
        }
        pi1_pi[k] = calculatePdfA(vs, k - 1, k) / PA_x0;
        debug_p[k] = std::make_pair(calculatePdfA(vs, k - 1, k), PA_x0);

        std::vector<double> p(t + s + 1);

        p[s] = 1;  // pdf for (s,t) path
        // this is not necessary
        //         if (path.eyeSize() > 0)
        //         {
        //             p[s] *= path.eyeEnd().getTotalPdfA();
        //         }
        //         if (path.lightSize() > 0)
        //         {
        //             p[s] *= path.lightEnd().getTotalPdfA();
        //         }

        for (int i = s; i <= k; i++)
        {
            p[i + 1] = p[i] * pi1_pi[i];
        }
        for (int i = s - 1; i >= 0; i--)
        {
            p[i] = p[i + 1] / pi1_pi[i];
        }

        for (int i = 0; i < vs.size() - 1; i++)
        {
            if (todo_getMaterial(vs[i]->obj->materialID()).isDelta())
            {
                p[i] = 0;
                p[i + 1] = 0;
            }
        }

        double mis_weight = 0;
        for (int i = 0; i < p.size(); i++)
        {
            double v = p[i] / p[s];
            mis_weight += v * v;
        }

        if (mis_weight != mis_weight)
        {
            return 0;
        }

        //         if (s == 0 && t == 6 && mis_weight < 10)
        //         {
        //             cout << "[ " << s << ", " << t << " ]  " << mis_weight <<
        //             endl;
        //         }

        //         if (s == 0)
        //             return 0;

        return 1.0 / mis_weight;
    }

    //////////////////////////////////////////////////////////////////////////

    vec3d pathThroughputWeighted(const CombinedPath& path)
    {
        vec3 f = vec3(1.0, 1.0, 1.0);

        // eye path inner weight
        for (int i = 1; i <= path.eyeSize() - 2; i++)
        {
            const Material& mat =
                todo_getMaterial(path.eye(i).obj->materialID());
            const vec3 d0 = normalize(path.eye(i - 1).p - path.eye(i).p);
            const vec3 d1 = normalize(path.eye(i + 1).p - path.eye(i).p);
            vec3       bsdf_weight = vec3(0.0);  // reflecting brdf
            if (!mat.isLight())
            {
                bsdf_weight = path.eye(i).weight;  // includes color
                bsdf_weight *= path.eye(i).cos_wo;
                bsdf_weight *= 1.0f / path.eye(i).pWo;
            }
            f *= bsdf_weight;
        }

        // light path inner weight
        for (int j = 1; j <= path.lightSize() - 2; j++)
        {
            const Material& mat =
                todo_getMaterial(path.light(j).obj->materialID());
            const vec3 d0 = normalize(path.light(j - 1).p - path.light(j).p);
            const vec3 d1 = normalize(path.light(j + 1).p - path.light(j).p);
            vec3       bsdf_weight = vec3(0.0);  // reflecting brdf
            if (!mat.isLight())
            {
                bsdf_weight = path.light(j).weight;  // includes color
                bsdf_weight *= path.light(j).cos_wo;
                bsdf_weight *= 1.0f / path.light(j).pWo;
            }
            f *= bsdf_weight;
        }

        //////////////////////////////////////////////////////////////////////////
        if (path.eyeSize() == 1)
        {
            f *= camera.measurementFunction(path.lightEnd());
        }

        if (path.lightSize() == 0)
        {
            if (path.eyeSize() > 1)
            {
                const Material& mat_eye_end =
                    todo_getMaterial(path.eyeEnd().obj->materialID());
                vec3 L = vec3(0);
                if (mat_eye_end.isLight())
                {
                    vec3 dir = path.eye(path.eyeSize() - 2).p - path.eyeEnd().p;
                    if (dot(dir, path.eyeEnd().n) > 0)
                    {
                        //                     L = mat_eye_end.color() / NUM_PI;
                        //                     // light
                        L = mat_eye_end.color();  // light
                    }
                }
                f *= L;
            }
        }

        if (path.lightSize() == 1)
        {
            float pdfA = 1.0f / scene.lightArea();
            f *= 1.0f / pdfA;

            if (path.eyeSize() > 1)
            {
                const Material& mat_eye_end =
                    todo_getMaterial(path.eyeEnd().obj->materialID());
                vec3 bsdf_eye_end = vec3(0.0);  // reflecting brdf
                if (!mat_eye_end.isLight())
                {
                    const vec3 d0 = normalize(path.eye(path.eyeSize() - 2).p -
                                              path.eyeEnd().p);
                    const vec3 d1 =
                        normalize(path.lightEnd().p - path.eyeEnd().p);
                    bsdf_eye_end = mat_eye_end.evaluate(
                        d0, path.eyeEnd().n, d1);  // includes color
                    bsdf_eye_end *=
                        geometryTerm(path.eyeEnd(), path.lightEnd());
                }
                f *= bsdf_eye_end;

                const Material& mat_light_end =
                    todo_getMaterial(path.lightEnd().obj->materialID());
                vec3 L = vec3(0);
                if (mat_light_end.isLight())
                {
                    vec3 dir = path.light(1).p - path.light(0).p;
                    if (dot(dir, path.light(0).n) > 0)
                    {
                        //                     L = mat_light_end.color() /
                        //                     NUM_PI; // light
                        L = mat_light_end.color();  // light
                    }
                }
                f *= L;
            }
        }

        if (path.lightSize() > 1)
        {
            const Material& mat_light_end =
                todo_getMaterial(path.lightEnd().obj->materialID());
            vec3 bsdf_light_end = vec3(0.0);  // reflecting brdf
            if (!mat_light_end.isLight())
            {
                const vec3 d0 = normalize(path.light(path.lightSize() - 2).p -
                                          path.lightEnd().p);
                const vec3 d1 = normalize(path.eyeEnd().p - path.lightEnd().p);
                bsdf_light_end =
                    mat_light_end.evaluate(d0, path.lightEnd().n, d1);
            }
            f *= bsdf_light_end;

            const Material& mat_light =
                todo_getMaterial(path.light(0).obj->materialID());
            assert(mat_light.type() == LGHT);
            vec3 L = vec3(0);
            if (mat_light.isLight())
            {
                vec3 dir = path.light(1).p - path.light(0).p;
                if (dot(dir, path.light(0).n) > 0)
                {
                    //                 L = mat_light.color(); // light
                    L = mat_light.color() * NUM_PI;  // light (pdfW = cos / pi,
                                                     // pdfA = cos' / r^2 * cos
                                                     // / pi, thus cancelled
                                                     // with G, becomes 1 / pi)
                }
            }
            f *= L;
            float pdfA = 1.0f / scene.lightArea();
            f *= 1.0f / pdfA;

            if (path.eyeSize() > 1)
            {
                const Material& mat_eye_end =
                    todo_getMaterial(path.eyeEnd().obj->materialID());
                vec3 bsdf_eye_end = vec3(0.0);  // reflecting brdf
                if (!mat_eye_end.isLight())
                {
                    const vec3 d0 = normalize(path.eye(path.eyeSize() - 2).p -
                                              path.eyeEnd().p);
                    const vec3 d1 =
                        normalize(path.lightEnd().p - path.eyeEnd().p);

                    bsdf_eye_end = mat_eye_end.evaluate(
                        d0, path.eyeEnd().n, d1);  // includes color
                    bsdf_eye_end *=
                        geometryTerm(path.eyeEnd(), path.lightEnd());
                }
                f *= bsdf_eye_end;
            }
        }

        return vec3d(f);
    }

public:
    ReferenceBDPT(const Camera& camera, FrameBuffer& film, const Scene& scene)
        : camera(camera), film(film), scene(scene)
    {
        pyramid.resize(MaxEvents * MaxEvents);
        for (auto& m : pyramid)
        {
            m.resize(film.width(), film.height());
        }
    }

    void render(const int MaxEyeEvents, StandardSampler& sampler)
    {
        // sample the path
        Path             eye_path = GenerateEyePath(MaxEyeEvents, sampler);
        Path             light_path = GenerateLightPath(MaxEyeEvents, sampler);
        PathContribution pc = CombinePaths(eye_path, light_path);

        // accumulate samples
        if (pc.sc > 0.0)
        {
            AccumulatePathContribution(pc);
        }
    }

    void dumpPyramid(float scale)
    {
        system("mkdir pyramid");
        for (int s = 0; s < MaxEvents; s++)
        {
            for (int t = 1; t < MaxEvents; t++)
            {
                if (s + t > MaxEvents)
                {
                    continue;
                }

                auto& img = pyramid[s * MaxEvents + t];

                std::string filename = std::string("pyramid/bdpt_s") + str(s) +
                                       "_t" + str(t) + ".ppm";
                img.scale(scale);
                img.tonemapGamma(2.2f);
                img.dumpPPM(filename.c_str());
            }
        }
    }
};

#endif  // bdpt_h__
