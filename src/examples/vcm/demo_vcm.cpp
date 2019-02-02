// Disclaimer: This demo is adapted from smallvcm http://www.smallvcm.com/

#include <omp.h>
#include <vector>
#include <unordered_map>
#include <memory>
#include <iostream>
#include <algorithm>
#include "vcm_material.h"
#include "vcm_scene.h"
#include "vcm_misc.h"
#include "vectors.h"
#include "intersection.h"
#include "camera.h"
#include "numeric.h"
#include "film.h"
#include "consoledebug.h"
#include "light.h"
#include "constants.h"

class DemoVCM
{
public:
    void run();
};

class HashGrid
{
public:
    void Reserve(int _numCells)
    {
        m_cellEnds.resize(_numCells);
    }

    template<typename tParticle>
    void Build(
        const std::vector<tParticle> &_particles,
        float _radius)
    {
        m_radius = _radius;
        m_radiusSqr = sq(m_radius);
        m_cellSize = m_radius * 2.0f;
        m_invCellSize = 1.0f / m_cellSize;

        m_BBoxMin = vec3(NUM_INFINITY);
        m_BBoxMax = vec3(-NUM_INFINITY);

        for (size_t i = 0; i < _particles.size(); i++)
        {
            const vec3 &pos = _particles[i].GetPosition();
            m_BBoxMax.x = f_max(m_BBoxMax.x, pos.x);
            m_BBoxMin.x = f_min(m_BBoxMin.x, pos.x);
            m_BBoxMax.y = f_max(m_BBoxMax.y, pos.y);
            m_BBoxMin.y = f_min(m_BBoxMin.y, pos.y);
            m_BBoxMax.z = f_max(m_BBoxMax.z, pos.z);
            m_BBoxMin.z = f_min(m_BBoxMin.z, pos.z);
        }

        m_indices.resize(_particles.size());
        memset(&m_cellEnds[0], 0, m_cellEnds.size() * sizeof(int));

        // set mCellEnds[x] to number of particles within x
        for (size_t i = 0; i < _particles.size(); i++)
        {
            const vec3 &pos = _particles[i].GetPosition();
            m_cellEnds[GetCellIndex(pos)]++;
        }

        // run exclusive prefix sum to really get the cell starts
        // mCellEnds[x] is now where the cell starts
        int sum = 0;
        for (size_t i = 0; i < m_cellEnds.size(); i++)
        {
            int temp = m_cellEnds[i];
            m_cellEnds[i] = sum;
            sum += temp;
        }

        for (size_t i = 0; i < _particles.size(); i++)
        {
            const vec3 &pos = _particles[i].GetPosition();
            const int targetIdx = m_cellEnds[GetCellIndex(pos)]++;
            m_indices[targetIdx] = int(i);
        }

        // now mCellEnds[x] points to the index right after the last
        // element of cell x

        //// DEBUG
        //for(size_t i=0; i<aParticles.size(); i++)
        //{
        //    const Vec3f &pos  = aParticles[i].GetPosition();
        //    Vec2i range = GetCellRange(GetCellIndex(pos));
        //    bool found = false;
        //    for(;range.x < range.y; range.x++)
        //    {
        //        if(mIndices[range.x] == i)
        //            found = true;
        //    }
        //    if(!found)
        //        printf("Error at particle %d\n", i);
        //}
    }

    template<typename tParticle, typename tQuery>
    void Process(
        const std::vector<tParticle> &_particles,
        tQuery& _query)
    {
        const vec3 queryPos = _query.GetPosition();

        const vec3 distMin = queryPos - m_BBoxMin;
        const vec3 distMax = m_BBoxMax - queryPos;
        if (distMin.x < 0.f) return;
        if (distMax.x < 0.f) return;
        if (distMin.y < 0.f) return;
        if (distMax.y < 0.f) return;
        if (distMin.z < 0.f) return;
        if (distMax.z < 0.f) return;

        const vec3 cellPt = m_invCellSize * distMin;
        const vec3 coordF(
            std::floor(cellPt.x),
            std::floor(cellPt.y),
            std::floor(cellPt.z));

        const int  px = int(coordF.x);
        const int  py = int(coordF.y);
        const int  pz = int(coordF.z);

        const vec3 fractCoord = cellPt - coordF;

        // cell size is twice the radius to ensure coverage
        const int  pxo = px + (fractCoord.x < 0.5f ? -1 : +1);
        const int  pyo = py + (fractCoord.y < 0.5f ? -1 : +1);
        const int  pzo = pz + (fractCoord.z < 0.5f ? -1 : +1);

        int found = 0;

        for (int j = 0; j < 8; j++)
        {
            int2 activeRange;
            switch (j)
            {
            case 0: activeRange = GetCellRange(GetCellIndex(int3(px, py, pz))); break;
            case 1: activeRange = GetCellRange(GetCellIndex(int3(px, py, pzo))); break;
            case 2: activeRange = GetCellRange(GetCellIndex(int3(px, pyo, pz))); break;
            case 3: activeRange = GetCellRange(GetCellIndex(int3(px, pyo, pzo))); break;
            case 4: activeRange = GetCellRange(GetCellIndex(int3(pxo, py, pz))); break;
            case 5: activeRange = GetCellRange(GetCellIndex(int3(pxo, py, pzo))); break;
            case 6: activeRange = GetCellRange(GetCellIndex(int3(pxo, pyo, pz))); break;
            case 7: activeRange = GetCellRange(GetCellIndex(int3(pxo, pyo, pzo))); break;
            }

            for (; activeRange.x < activeRange.y; activeRange.x++)
            {
                const int particleIndex = m_indices[activeRange.x];
                const tParticle &particle = _particles[particleIndex];

                const float distSqr =
                    lengthSquared(_query.GetPosition() - particle.GetPosition());

                if (distSqr <= m_radiusSqr)
                    _query.Process(particle);
            }
        }
    }

private:

    int2 GetCellRange(int _cellIndex) const
    {
        if (_cellIndex == 0)
        {
            return int2(0, m_cellEnds[0], 0);
        }
        return int2(m_cellEnds[_cellIndex - 1], m_cellEnds[_cellIndex], 0);
    }

    int GetCellIndex(const int3 &_coord) const
    {
        uint32_t x = uint32_t(_coord.x);
        uint32_t y = uint32_t(_coord.y);
        uint32_t z = uint32_t(_coord.z);

        return int(((x * 73856093) ^ (y * 19349663) ^
            (z * 83492791)) % uint32_t(m_cellEnds.size()));
    }

    int GetCellIndex(const vec3 &_point) const
    {
        const vec3 distMin = _point - m_BBoxMin;

        const vec3 coordF(
            std::floor(m_invCellSize * distMin.x),
            std::floor(m_invCellSize * distMin.y),
            std::floor(m_invCellSize * distMin.z));

        const int3 coordI = int3(int(coordF.x), int(coordF.y), int(coordF.z));

        return GetCellIndex(coordI);
    }

private:

    vec3 m_BBoxMin;
    vec3 m_BBoxMax;
    std::vector<int> m_indices;
    std::vector<int> m_cellEnds;

    float m_radius;
    float m_radiusSqr;
    float m_cellSize;
    float m_invCellSize;
};

////////////////////////////////////////////////////////////////////////////////
// A NOTE ON PATH MIS WEIGHT EVALUATION
////////////////////////////////////////////////////////////////////////////////
//
// We compute path MIS weights iteratively as we trace the light and eye
// sub-paths. We cache three floating points quantities at each sub-path vertex:
//
//   dVCM  dVC  dVM
//
// These quantities represent partial weights associated with the sub-path. When
// we connect or merge one vertex to another, we use these quantities to quickly
// evaluate the MIS weight for the full path we have constructed. This scheme is
// presented in the technical report
//
//   "Implementing Vertex Connection and Merging"
//   http://www.iliyan.com/publications/ImplementingVCM
//
// The MIS code in the VertexCM class references the corresponding equations in
// the report in the form
//
//   [tech. rep. (##)]
//
// where ## is the equation number. 
//

class Pyramid
{
    std::unordered_map<int, std::unique_ptr<FrameBuffer>> _pyramid;
    int _w, _h;
    int _max_length;

    static int hash(int s, int y)
    {
        return ((s & 0xffff) << 16) + (y & 0xffff);
    }
    static void unhash(int id, int& s, int& t)
    {
        s = (id >> 16) & 0xffff;
        t = id & 0xffff;
    }

public:
    Pyramid(int w, int h, int max_length)
        : _w(w), _h(h), _max_length(max_length)
    {
        system("mkdir pyramid2");
    }
    std::unique_ptr<FrameBuffer>& get(int s, int t)
    {
        int id = hash(s, t);
        if (_pyramid.find(id) == _pyramid.end())
        {
            _pyramid.insert(std::make_pair(id, std::unique_ptr<FrameBuffer>(new FrameBuffer(_w, _h))));
        }
        return _pyramid[id];
    }
    void dump(float scale)
    {
        for (auto& p : _pyramid)
        {
            int s, t;
            unhash(p.first, s, t);
            if (s + t > _max_length)
            {
                continue;
            }

            std::string filename = std::string("pyramid2/bdpt_s") + str(s) + "_t" + str(t) + ".ppm";
            p.second->scale(scale);
            p.second->tonemapGamma(2.2f);
            p.second->dumpPPM(filename.c_str());
        }
    }
};

class VertexCM
{
public:
    uint32_t            m_maxPathLength;
    uint32_t            m_minPathLength;

    uint32_t            m_iterations;
    FrameBuffer         m_framebuffer;
    const Scene_debug&  m_scene;

    std::unique_ptr<Pyramid> m_pyramid;

    // final output of this rendering work thread
    void GetFramebuffer(FrameBuffer& framebuffer_)
    {
        framebuffer_ = m_framebuffer;

        if (m_iterations > 0)
        {
            framebuffer_.scale(1.0f / m_iterations);
        }
    }

    void dumpPyramid()
    {
        if (m_iterations > 0)
        {
            m_pyramid->dump(1.0f / m_iterations);
        }
    }

private:
    // The sole point of this structure is to make carrying around the ray baggage easier.
    struct SubPathState
    {
        vec3 mOrigin;             // Path origin
        vec3 mDirection;          // Where to go next
        vec3 mThroughput;         // Path throughput
        uint32_t mPathLength : 30; // Number of path segments, including this
        uint32_t mIsFiniteLight : 1; // Just generate by finite light
        uint32_t mSpecularPath : 1; // All scattering events so far were specular

        float dVCM; // MIS quantity used for vertex connection and merging
        float dVC;  // MIS quantity used for vertex connection
        float dVM;  // MIS quantity used for vertex merging
    };

    // Path vertex, used for merging and connection
    template<bool tFromLight>
    struct PathVertex
    {
        vec3 mHitpoint;   // Position of the vertex
        vec3 mThroughput; // Path throughput (including emission)
        uint32_t  mPathLength; // Number of segments between source and vertex

        // Stores all required local information, including incoming direction.
        BSDF<tFromLight> mBsdf;

        float dVCM; // MIS quantity used for vertex connection and merging
        float dVC;  // MIS quantity used for vertex connection
        float dVM;  // MIS quantity used for vertex merging

        // Used by HashGrid
        const vec3& GetPosition() const
        {
            return mHitpoint;
        }
    };

    typedef PathVertex<false> CameraVertex;
    typedef PathVertex<true>  LightVertex;

    typedef BSDF<false>       CameraBSDF;
    typedef BSDF<true>        LightBSDF;

    // Range query used for PPM, BPT, and VCM. When HashGrid finds a vertex
    // within range -- Process() is called and vertex
    // merging is performed. BSDF of the camera vertex is used.
    class RangeQuery
    {
    public:

        RangeQuery(
            const VertexCM     &aVertexCM,
            const vec3        &aCameraPosition,
            const CameraBSDF   &aCameraBsdf,
            const SubPathState &aCameraState
            ) :
            mVertexCM(aVertexCM),
            mCameraPosition(aCameraPosition),
            mCameraBsdf(aCameraBsdf),
            mCameraState(aCameraState),
            mContrib(0)
        {}

        const vec3& GetPosition() const { return mCameraPosition; }

        const vec3& GetContrib() const { return mContrib; }

        void Process(const LightVertex& aLightVertex)
        {
            // Reject if full path length below/above min/max path length
            if ((aLightVertex.mPathLength + mCameraState.mPathLength > mVertexCM.m_maxPathLength) ||
                (aLightVertex.mPathLength + mCameraState.mPathLength < mVertexCM.m_minPathLength))
            {
                return;
            }

            // Retrieve light incoming direction in world coordinates
            const vec3 lightDirection = aLightVertex.mBsdf.WorldDirFix();

            float cosCamera, cameraBsdfDirPdfW, cameraBsdfRevPdfW;
            const vec3 cameraBsdfFactor = mCameraBsdf.Evaluate(
                lightDirection, cosCamera, &cameraBsdfDirPdfW,
                &cameraBsdfRevPdfW);

            if (isZero(cameraBsdfFactor))
            {
                return;
            }

            cameraBsdfDirPdfW *= mCameraBsdf.ContinuationProb();

            // Even though this is pdf from camera BSDF, the continuation probability
            // must come from light BSDF, because that would govern it if light path
            // actually continued
            cameraBsdfRevPdfW *= aLightVertex.mBsdf.ContinuationProb();

            // Partial light sub-path MIS weight [tech. rep. (38)]
            const float wLight = aLightVertex.dVCM * mVertexCM.mMisVcWeightFactor +
                aLightVertex.dVM * mVertexCM.Mis(cameraBsdfDirPdfW);

            // Partial eye sub-path MIS weight [tech. rep. (39)]
            const float wCamera = mCameraState.dVCM * mVertexCM.mMisVcWeightFactor +
                mCameraState.dVM * mVertexCM.Mis(cameraBsdfRevPdfW);

            // Full path MIS weight [tech. rep. (37)]. No MIS for PPM
            const float misWeight = mVertexCM.mPpm ?
                1.f :
                1.f / (wLight + 1.f + wCamera);

            mContrib += misWeight * cameraBsdfFactor * aLightVertex.mThroughput;
        }

    private:

        const VertexCM&         mVertexCM;
        const vec3&             mCameraPosition;
        const CameraBSDF&       mCameraBsdf;
        const SubPathState&     mCameraState;
        vec3                    mContrib;
    };

private:
    bool  mUseVM;             // Vertex merging (of some form) is used
    bool  mUseVC;             // Vertex connection (BPT) is used
    bool  mLightTraceOnly;    // Do only light tracing
    bool  mPpm;               // Do PPM, same terminates camera after first merge

    float mRadiusAlpha;       // Radius reduction rate parameter
    float mBaseRadius;        // Initial merging radius
    float mMisVmWeightFactor; // Weight of vertex merging (used in VC)
    float mMisVcWeightFactor; // Weight of vertex connection (used in VM)
    float mScreenPixelCount;  // Number of pixels
    float mLightSubPathCount; // Number of light sub-paths
    float mVmNormalization;   // 1 / (Pi * radius^2 * light_path_count)

    std::vector<LightVertex> mLightVertices; //!< Stored light vertices

    // For light path belonging to pixel index [x] it stores
    // where it's light vertices end (begin is at [x-1])
    std::vector<int> mPathEnds;
    HashGrid         mHashGrid;

public:
    enum AlgorithmType
    {
        // light vertices contribute to camera,
        // No MIS weights (dVCM, dVM, dVC all ignored)
        kLightTrace = 0,

        // Camera and light vertices merged on first non-specular surface from camera.
        // Cannot handle mixed specular + non-specular materials.
        // No MIS weights (dVCM, dVM, dVC all ignored)
        kPpm,

        // Camera and light vertices merged on along full path.
        // dVCM and dVM used for MIS
        kBpm,

        // Standard bidirectional path tracing
        // dVCM and dVC used for MIS
        kBpt,

        // Vertex connection and mering
        // dVCM, dVM, and dVC used for MIS
        kVcm
    };

public:
    VertexCM(
        const Scene_debug&  aScene,
        AlgorithmType aAlgorithm,
        const float   aRadiusFactor,
        const float   aRadiusAlpha
        ) :
        m_scene(aScene),
        mLightTraceOnly(false),
        mUseVC(false),
        mUseVM(false),
        mPpm(false)
    {
        m_minPathLength = 0;
        m_maxPathLength = 10; // 20;// 2;
        m_iterations = 0;
        m_framebuffer.resize(m_scene.m_camera->filmWidth(), m_scene.m_camera->filmHeight());

        switch (aAlgorithm)
        {
        case kLightTrace:
            mLightTraceOnly = true;
            break;
        case kPpm:
            mPpm = true;
            mUseVM = true;
            break;
        case kBpm:
            mUseVM = true;
            break;
        case kBpt:
            mUseVC = true;
            break;
        case kVcm:
            mUseVC = true;
            mUseVM = true;
            break;
        default:
            printf("Unknown algorithm requested\n");
            break;
        }

        if (mPpm)
        {
            // We will check the scene to make sure it does not contain mixed
            // specular and non-specular materials
            for (int i = 0; i < m_scene.GetMaterialCount(); ++i)
            {
                const Material_debug &mat = m_scene.GetMaterial(i);

                const bool hasNonSpecular =
                    (maxComponent(mat.m_diffuseReflectance) > 0) ||
                    (maxComponent(mat.m_phongReflectance) > 0);

                const bool hasSpecular =
                    (maxComponent(mat.m_mirrorReflectance) > 0) ||
                    (mat.m_IOR > 0);

                if (hasNonSpecular && hasSpecular)
                {
                    printf(
                        "*WARNING* Our PPM implementation cannot handle materials mixing\n"
                        "Specular and NonSpecular BSDFs. The extension would be\n"
                        "fairly straightforward. In SampleScattering for camera sub-paths\n"
                        "limit the considered events to Specular only.\n"
                        "Merging will use non-specular components, scattering will be specular.\n"
                        "If there is no specular component, the ray will terminate.\n\n");

                    printf("We are now switching from *PPM* to *BPM*, which can handle the scene\n\n");

                    mPpm = false;
                    break;
                }
            }
        }

        mBaseRadius = aRadiusFactor * m_scene.m_sceneSphere.m_radius;
        std::cout << "!!base radius = " << mBaseRadius << std::endl;

        mRadiusAlpha = aRadiusAlpha;

        m_pyramid = std::make_unique<Pyramid>(m_scene.m_camera->filmWidth(), m_scene.m_camera->filmHeight(), m_maxPathLength);
    }

    virtual void RunIteration(int aIteration)
    {
        // While we have the same number of pixels (camera paths)
        // and light paths, we do keep them separate for clarity reasons
        const int resX = int(m_scene.m_camera->filmWidth());
        const int resY = int(m_scene.m_camera->filmHeight());
        const int pathCount = resX * resY;
        mScreenPixelCount = float(resX * resY);
        mLightSubPathCount = float(resX * resY);

        // Setup our radius, 1st iteration has aIteration == 0, thus offset
        float radius = mBaseRadius;
        radius /= std::pow(float(aIteration + 1), 0.5f * (1 - mRadiusAlpha));
        // Purely for numeric stability
        radius = f_max(radius, 1e-7f);
        const float radiusSqr = sq(radius);

        // Factor used to normalise vertex merging contribution.
        // We divide the summed up energy by disk radius and number of light paths
        mVmNormalization = 1.0f / (radiusSqr * NUM_PI * mLightSubPathCount);

        // MIS weight constant [tech. rep. (20)], with n_VC = 1 and n_VM = mLightPathCount
        const float etaVCM = (NUM_PI * radiusSqr) * mLightSubPathCount;
        mMisVmWeightFactor = mUseVM ? Mis(etaVCM) : 0.f;
        mMisVcWeightFactor = mUseVC ? Mis(1.f / etaVCM) : 0.f;

        // Clear path ends, nothing ends anywhere
        mPathEnds.resize(pathCount);
        memset(&mPathEnds[0], 0, mPathEnds.size() * sizeof(int));

        // Remove all light vertices and reserve space for some
        mLightVertices.reserve(pathCount);
        mLightVertices.clear();

        //////////////////////////////////////////////////////////////////////////
        // Generate light paths
        //////////////////////////////////////////////////////////////////////////
        for (int pathIdx = 0; pathIdx < pathCount; pathIdx++)
        {
            SubPathState lightState;
            GenerateLightSample(lightState);

            //////////////////////////////////////////////////////////////////////////
            // Trace light path
            for (;; ++lightState.mPathLength)
            {
                // Offset ray origin instead of setting tmin due to numeric
                // issues in ray-sphere intersection. The isect.dist has to be
                // extended by this EPS_RAY after hit point is determined
                Ray ray(lightState.mOrigin /*+ lightState.mDirection * NUM_EPS_RAY*/,
                    lightState.mDirection);

                HitInfo_debug hit = m_scene.intersect(ray);
                if (!hit)
                {
                    break;
                }

                const vec3 hitPoint = ray.at(hit.distance());
                const Material_debug& mat = m_scene.m_materials[hit.materialID()];

                LightBSDF bsdf(ray, hit.shadingNormal(), &mat);
                if (!bsdf.IsValid())
                {
                    break;
                }

                // Update the MIS quantities before storing them at the vertex.
                // These updates follow the initialization in GenerateLightSample() or
                // SampleScattering(), and together implement equations [tech. rep. (31)-(33)]
                // or [tech. rep. (34)-(36)], respectively.
                {
                    // Infinite lights use MIS handled via solid angle integration,
                    // so do not divide by the distance for such lights [tech. rep. Section 5.1]
                    if (lightState.mPathLength > 1 || lightState.mIsFiniteLight == 1)
                    {
                        lightState.dVCM *= Mis(sq(hit.distance()));
                    }

                    lightState.dVCM /= Mis(std::abs(bsdf.CosThetaFix()));
                    lightState.dVC /= Mis(std::abs(bsdf.CosThetaFix()));
                    lightState.dVM /= Mis(std::abs(bsdf.CosThetaFix()));
                }

                // Store vertex, unless BSDF is purely specular, which prevents
                // vertex connections and merging
                if (!bsdf.IsDelta() && (mUseVC || mUseVM))
                {
                    LightVertex lightVertex;
                    lightVertex.mHitpoint = hitPoint;
                    lightVertex.mThroughput = lightState.mThroughput;
                    lightVertex.mPathLength = lightState.mPathLength;
                    lightVertex.mBsdf = bsdf;

                    lightVertex.dVCM = lightState.dVCM;
                    lightVertex.dVC = lightState.dVC;
                    lightVertex.dVM = lightState.dVM;

                    mLightVertices.push_back(lightVertex);
                }

                // Connect to camera, unless BSDF is purely specular
                if (!bsdf.IsDelta() && (mUseVC || mLightTraceOnly))
                {
                    if (lightState.mPathLength + 1 >= m_minPathLength)
                    {
                        ConnectToCamera(lightState, hitPoint, bsdf);
                    }
                }

                // Terminate if the path would become too long after scattering
                if (lightState.mPathLength + 2 > m_maxPathLength)
                {
                    break;
                }

                // Continue random walk
                if (!SampleScattering(bsdf, hitPoint, lightState))
                {
                    break;
                }

                if (pathIdx == 10000)
                {
                    std::cout << "[ s = " << lightState.mPathLength << " ]" << std::endl;
                }
            }

            mPathEnds[pathIdx] = (int)mLightVertices.size();
        }

        //////////////////////////////////////////////////////////////////////////
        // Build hash grid
        //////////////////////////////////////////////////////////////////////////

        // Only build grid when merging (VCM, BPM, and PPM)
        if (mUseVM)
        {
            // The number of cells is somewhat arbitrary, but seems to work ok
            mHashGrid.Reserve(pathCount);
            mHashGrid.Build(mLightVertices, radius);

            std::cout << "!!light verts " << mLightVertices.size() << std::endl;
        }

        //////////////////////////////////////////////////////////////////////////
        // Generate camera paths
        //////////////////////////////////////////////////////////////////////////

        // Unless rendering with traditional light tracing
        for (int pathIdx = 0; (pathIdx < pathCount) && (!mLightTraceOnly); ++pathIdx)
        {
            SubPathState cameraState;
            const vec2 screenSample = GenerateCameraSample(pathIdx, cameraState);
            vec3 color(0);

            //////////////////////////////////////////////////////////////////////
            // Trace camera path
            for (;; ++cameraState.mPathLength)
            {
                // Offset ray origin instead of setting tmin due to numeric
                // issues in ray-sphere intersection. The isect.dist has to be
                // extended by this EPS_RAY after hit point is determined
                Ray ray(cameraState.mOrigin/* + cameraState.mDirection * NUM_EPS_RAY*/,
                    cameraState.mDirection);

                HitInfo_debug hit = m_scene.intersect(ray);

                // Get radiance from environment
                if (!hit)
                {
                    if (m_scene.GetBackground() != nullptr)
                    {
                        if (cameraState.mPathLength >= m_minPathLength)
                        {
                            vec3 _color = cameraState.mThroughput *
                                GetLightRadiance(m_scene.GetBackground(), cameraState,
                                vec3(0), ray.dir);
                            color += _color;
                            if (m_pyramid)
                            {
                                m_pyramid->get(0, cameraState.mPathLength + 1)->accumulatePixel((int)screenSample.x, (int)screenSample.y, _color);
                            }
                        }
                    }

                    break;
                }

                const vec3 hitPoint = ray.at(hit.distance());
                const Material_debug& mat = m_scene.m_materials[hit.materialID()];

                CameraBSDF bsdf(ray, hit.shadingNormal(), &mat);
                if (!bsdf.IsValid())
                {
                    break;
                }

                // Update the MIS quantities, following the initialization in
                // GenerateLightSample() or SampleScattering(). Implement equations
                // [tech. rep. (31)-(33)] or [tech. rep. (34)-(36)], respectively.
                {
                    cameraState.dVCM *= Mis(sq(hit.distance()));
                    cameraState.dVCM /= Mis(std::abs(bsdf.CosThetaFix()));
                    cameraState.dVC /= Mis(std::abs(bsdf.CosThetaFix()));
                    cameraState.dVM /= Mis(std::abs(bsdf.CosThetaFix()));
                }

                int light_id = m_scene.toLightID(hit.materialID());
                // Light source has been hit; terminate afterwards, since
                // our light sources do not have reflective properties
                if (light_id >= 0)
                {
                    const AbstractLight *light = m_scene.GetLightPtr(light_id);

                    if (cameraState.mPathLength >= m_minPathLength)
                    {
                        vec3 _color = cameraState.mThroughput *
                            GetLightRadiance(light, cameraState, hitPoint, ray.dir);
                        color += _color;
                        if (m_pyramid)
                        {
                            m_pyramid->get(0, cameraState.mPathLength + 1)->accumulatePixel((int)screenSample.x, (int)screenSample.y, _color);
                        }
                    }

                    break;
                }

                // Terminate if eye sub-path is too long for connections or merging
                if (cameraState.mPathLength >= m_maxPathLength)
                {
                    break;
                }

                ////////////////////////////////////////////////////////////////
                // Vertex connection: Connect to a light source
                if (!bsdf.IsDelta() && mUseVC)
                {
                    if (cameraState.mPathLength + 1 >= m_minPathLength)
                    {
                        vec3 _color = cameraState.mThroughput *
                            DirectIllumination(cameraState, hitPoint, bsdf);
                        color += _color;
                        if (m_pyramid)
                        {
                            m_pyramid->get(1, cameraState.mPathLength + 1)->accumulatePixel((int)screenSample.x, (int)screenSample.y, _color);
                        }
                    }
                }

                ////////////////////////////////////////////////////////////////
                // Vertex connection: Connect to light vertices
                if (!bsdf.IsDelta() && mUseVC)
                {
                    // For VC, each light sub-path is assigned to a particular eye
                    // sub-path, as in traditional BPT. It is also possible to
                    // connect to vertices from any light path, but MIS should
                    // be revisited.
                    const int2 range(
                        (pathIdx == 0) ? 0 : mPathEnds[pathIdx - 1],
                        mPathEnds[pathIdx], 0);

                    for (int i = range.x; i < range.y; i++)
                    {
                        const LightVertex &lightVertex = mLightVertices[i];

                        if (lightVertex.mPathLength + 1 +
                            cameraState.mPathLength < m_minPathLength)
                        {
                            continue;
                        }

                        // Light vertices are stored in increasing path length
                        // order; once we go above the max path length, we can
                        // skip the rest
                        if (lightVertex.mPathLength + 1 +
                            cameraState.mPathLength > m_maxPathLength)
                        {
                            break;
                        }

                        vec3 _color = cameraState.mThroughput * lightVertex.mThroughput *
                            ConnectVertices(lightVertex, bsdf, hitPoint, cameraState);
                        color += _color;
                        if (m_pyramid)
                        {
                            m_pyramid->get(lightVertex.mPathLength + 1, cameraState.mPathLength + 1)->accumulatePixel((int)screenSample.x, (int)screenSample.y, _color);
                        }
                    }
                }

                ////////////////////////////////////////////////////////////////
                // Vertex merging: Merge with light vertices
                if (!bsdf.IsDelta() && mUseVM)
                {
//                     std::cout << "!!merging" << std::endl;

                    RangeQuery query(*this, hitPoint, bsdf, cameraState);
                    mHashGrid.Process(mLightVertices, query);
                    vec3 _color = cameraState.mThroughput * mVmNormalization * query.GetContrib();
                    color += _color;

                    // PPM merges only at the first non-specular surface from camera
                    if (mPpm)
                    {
                        break;
                    }
                }

                if (!SampleScattering(bsdf, hitPoint, cameraState))
                {
                    break;
                }

                if (pathIdx == 10000)
                {
                    std::cout << "[ t = " << cameraState.mPathLength << " ]" << std::endl;
                }
            }

            m_framebuffer.accumulatePixel((int)screenSample.x, (int)screenSample.y, color);
        }

        m_iterations++;
    }

private:

    // Mis power, we use balance heuristic
    float Mis(float aPdf) const
    {
        //return std::pow(aPdf, /*power*/);
        return aPdf;
    }

    //////////////////////////////////////////////////////////////////////////
    // Camera tracing methods
    //////////////////////////////////////////////////////////////////////////

    // Generates new camera sample given a pixel index
    vec2 GenerateCameraSample(
        const int    aPixelIndex,
        SubPathState &oCameraState)
    {
        const Camera &camera = *m_scene.m_camera;
        const int resX = int(camera.filmWidth());
        const int resY = int(camera.filmHeight());

        // Determine pixel (x, y)
        const int x = aPixelIndex % resX;
        const int y = aPixelIndex / resY;

        // Jitter pixel position
        const vec2 sample = vec2(float(x), float(y), 0) + vec2(randf(), randf(), 0);

        // Generate ray
        const Ray primaryRay = camera.makeRay(sample.x / resX, sample.y / resY);

        // Compute pdf conversion factor from area on image plane to solid angle on ray
        const float cosAtCamera = dot(camera.forwardDirection(), primaryRay.dir);
        const float imagePointToCameraDist = camera.filmDistance() / cosAtCamera;
        const float imageToSolidAngleFactor = sq(imagePointToCameraDist) / cosAtCamera;

        // We put the virtual image plane at such a distance from the camera origin
        // that the pixel area is one and thus the image plane sampling pdf is 1.
        // The solid angle ray pdf is then equal to the conversion factor from
        // image plane area density to ray solid angle density
//         const float cameraPdfW = imageToSolidAngleFactor;
        const float cameraPdfW = camera.pixelPdfA() * imageToSolidAngleFactor;

        oCameraState.mOrigin = primaryRay.orig;
        oCameraState.mDirection = primaryRay.dir;
        oCameraState.mThroughput = vec3(1);

        oCameraState.mPathLength = 1;
        oCameraState.mSpecularPath = 1;

        // Eye sub-path MIS quantities. Implements [tech. rep. (31)-(33)] partially.
        // The evaluation is completed after tracing the camera ray in the eye sub-path loop.
        oCameraState.dVCM = Mis(mLightSubPathCount / cameraPdfW);
        oCameraState.dVC = 0;
        oCameraState.dVM = 0;

        return sample;
    }

    // Returns the radiance of a light source when hit by a random ray,
    // multiplied by MIS weight. Can be used for both Background and Area lights.
    //
    // For Background lights:
    //    Has to be called BEFORE updating the MIS quantities.
    //    Value of aHitpoint is irrelevant (passing Vec3f(0))
    //
    // For Area lights:
    //    Has to be called AFTER updating the MIS quantities.
    vec3 GetLightRadiance(
        const AbstractLight *aLight,
        const SubPathState  &aCameraState,
        const vec3         &aHitpoint,
        const vec3         &aRayDirection) const
    {
        // We sample lights uniformly
        const int   lightCount = m_scene.GetLightCount();
        const float lightPickProb = 1.f / lightCount;

        float directPdfA, emissionPdfW;
        const vec3 radiance = aLight->GetRadiance(m_scene.m_sceneSphere,
            aRayDirection, aHitpoint, &directPdfA, &emissionPdfW);

        if (isZero(radiance))
        {
            return vec3(0);
        }

        // If we see light source directly from camera, no weighting is required
        if (aCameraState.mPathLength == 1)
        {
            return radiance;
        }

        // When using only vertex merging, we want purely specular paths
        // to give radiance (cannot get it otherwise). Rest is handled
        // by merging and we should return 0.
        if (mUseVM && !mUseVC)
        {
            return aCameraState.mSpecularPath ? radiance : vec3(0);
        }

        directPdfA *= lightPickProb;
        emissionPdfW *= lightPickProb;

        // Partial eye sub-path MIS weight [tech. rep. (43)].
        // If the last hit was specular, then dVCM == 0.
        const float wCamera = Mis(directPdfA) * aCameraState.dVCM +
            Mis(emissionPdfW) * aCameraState.dVC;

        // Partial light sub-path weight is 0 [tech. rep. (42)].

        // Full path MIS weight [tech. rep. (37)].
        const float misWeight = 1.f / (1.f + wCamera);

        return misWeight * radiance;
    }

    // Connects camera vertex to randomly chosen light point.
    // Returns emitted radiance multiplied by path MIS weight.
    // Has to be called AFTER updating the MIS quantities.
    vec3 DirectIllumination(
        const SubPathState &aCameraState,
        const vec3        &aHitpoint,
        const CameraBSDF   &aBsdf)
    {
        // We sample lights uniformly
        const int   lightCount = m_scene.GetLightCount();
        const float lightPickProb = 1.f / lightCount;

        const int   lightID = int(randf() * lightCount);
        const vec2 rndPosSamples = vec2(randf(), randf(), 0);

        const AbstractLight *light = m_scene.GetLightPtr(lightID);

        vec3 directionToLight;
        float distance;
        float directPdfW, emissionPdfW, cosAtLight;
        const vec3 radiance = light->Illuminate(m_scene.m_sceneSphere, aHitpoint,
            rndPosSamples, directionToLight, distance, directPdfW,
            &emissionPdfW, &cosAtLight);

        // If radiance == 0, other values are undefined, so have to early exit
        if (isZero(radiance))
        {
            return vec3(0);
        }

        float bsdfDirPdfW, bsdfRevPdfW, cosToLight;
        const vec3 bsdfFactor = aBsdf.Evaluate(
            directionToLight, cosToLight, &bsdfDirPdfW, &bsdfRevPdfW);

        if (isZero(bsdfFactor))
        {
            return vec3(0);
        }

        const float continuationProbability = aBsdf.ContinuationProb();

        // If the light is delta light, we can never hit it
        // by BSDF sampling, so the probability of this path is 0
        bsdfDirPdfW *= light->IsDelta() ? 0.f : continuationProbability;

        bsdfRevPdfW *= continuationProbability;

        // Partial light sub-path MIS weight [tech. rep. (44)].
        // Note that wLight is a ratio of area pdfs. But since both are on the
        // light source, their distance^2 and cosine terms cancel out.
        // Therefore we can write wLight as a ratio of solid angle pdfs,
        // both expressed w.r.t. the same shading point.
        const float wLight = Mis(bsdfDirPdfW / (lightPickProb * directPdfW));

        // Partial eye sub-path MIS weight [tech. rep. (45)].
        //
        // In front of the sum in the parenthesis we have Mis(ratio), where
        //    ratio = emissionPdfA / directPdfA,
        // with emissionPdfA being the product of the pdfs for choosing the
        // point on the light source and sampling the outgoing direction.
        // What we are given by the light source instead are emissionPdfW
        // and directPdfW. Converting to area pdfs and plugging into ratio:
        //    emissionPdfA = emissionPdfW * cosToLight / dist^2
        //    directPdfA   = directPdfW * cosAtLight / dist^2
        //    ratio = (emissionPdfW * cosToLight / dist^2) / (directPdfW * cosAtLight / dist^2)
        //    ratio = (emissionPdfW * cosToLight) / (directPdfW * cosAtLight)
        //
        // Also note that both emissionPdfW and directPdfW should be
        // multiplied by lightPickProb, so it cancels out.
        const float wCamera = Mis(emissionPdfW * cosToLight / (directPdfW * cosAtLight)) * (
            mMisVmWeightFactor + aCameraState.dVCM + aCameraState.dVC * Mis(bsdfRevPdfW));

        // Full path MIS weight [tech. rep. (37)]
        const float misWeight = 1.f / (wLight + 1.f + wCamera);

        const vec3 contrib =
            (misWeight * cosToLight / (lightPickProb * directPdfW)) * (radiance * bsdfFactor);

        if (isZero(contrib) || m_scene.Occluded(aHitpoint, directionToLight, distance))
        {
            return vec3(0);
        }

        return contrib;
    }

    // Connects an eye and a light vertex. Result multiplied by MIS weight, but
    // not multiplied by vertex throughputs. Has to be called AFTER updating MIS
    // constants. 'direction' is FROM eye TO light vertex.
    vec3 ConnectVertices(
        const LightVertex  &aLightVertex,
        const CameraBSDF   &aCameraBsdf,
        const vec3        &aCameraHitpoint,
        const SubPathState &aCameraState) const
    {
        // Get the connection
        vec3 direction = aLightVertex.mHitpoint - aCameraHitpoint;
        const float dist2 = lengthSquared(direction);
        float  distance = std::sqrt(dist2);
        direction /= distance;

        // Evaluate BSDF at camera vertex
        float cosCamera, cameraBsdfDirPdfW, cameraBsdfRevPdfW;
        const vec3 cameraBsdfFactor = aCameraBsdf.Evaluate(
            direction, cosCamera, &cameraBsdfDirPdfW,
            &cameraBsdfRevPdfW);

        if (isZero(cameraBsdfFactor))
        {
            return vec3(0);
        }

        // Camera continuation probability (for Russian roulette)
        const float cameraCont = aCameraBsdf.ContinuationProb();
        cameraBsdfDirPdfW *= cameraCont;
        cameraBsdfRevPdfW *= cameraCont;

        // Evaluate BSDF at light vertex
        float cosLight, lightBsdfDirPdfW, lightBsdfRevPdfW;
        const vec3 lightBsdfFactor = aLightVertex.mBsdf.Evaluate(
            -direction, cosLight, &lightBsdfDirPdfW,
            &lightBsdfRevPdfW);

        if (isZero(lightBsdfFactor))
        {
            return vec3(0);
        }

        // Light continuation probability (for Russian roulette)
        const float lightCont = aLightVertex.mBsdf.ContinuationProb();
        lightBsdfDirPdfW *= lightCont;
        lightBsdfRevPdfW *= lightCont;

        // Compute geometry term
        const float geometryTerm = cosLight * cosCamera / dist2;
        if (geometryTerm < 0)
        {
            return vec3(0);
        }

        // Convert pdfs to area pdf
        const float cameraBsdfDirPdfA = pdfWtoA(cameraBsdfDirPdfW, distance, cosLight);
        const float lightBsdfDirPdfA = pdfWtoA(lightBsdfDirPdfW, distance, cosCamera);

        // Partial light sub-path MIS weight [tech. rep. (40)]
        const float wLight = Mis(cameraBsdfDirPdfA) * (
            mMisVmWeightFactor + aLightVertex.dVCM + aLightVertex.dVC * Mis(lightBsdfRevPdfW));

        // Partial eye sub-path MIS weight [tech. rep. (41)]
        const float wCamera = Mis(lightBsdfDirPdfA) * (
            mMisVmWeightFactor + aCameraState.dVCM + aCameraState.dVC * Mis(cameraBsdfRevPdfW));

        // Full path MIS weight [tech. rep. (37)]
        const float misWeight = 1.f / (wLight + 1.f + wCamera);

        const vec3 contrib = (misWeight * geometryTerm) * cameraBsdfFactor * lightBsdfFactor;

        if (isZero(contrib) || m_scene.Occluded(aCameraHitpoint, direction, distance))
        {
            return vec3(0);
        }

        return contrib;
    }

    //////////////////////////////////////////////////////////////////////////
    // Light tracing methods
    //////////////////////////////////////////////////////////////////////////

    // Samples light emission
    void GenerateLightSample(SubPathState &oLightState)
    {
        // We sample lights uniformly
        const int   lightCount = m_scene.GetLightCount();
        const float lightPickProb = 1.f / lightCount;

        const int   lightID = int(randf() * lightCount);
        const vec2 rndDirSamples = vec2(randf(), randf(), 0);
        const vec2 rndPosSamples = vec2(randf(), randf(), 0);

        const AbstractLight *light = m_scene.GetLightPtr(lightID);

        float emissionPdfW, directPdfW, cosLight;
        oLightState.mThroughput = light->Emit(m_scene.m_sceneSphere, rndDirSamples, rndPosSamples,
            oLightState.mOrigin, oLightState.mDirection,
            emissionPdfW, &directPdfW, &cosLight);

        emissionPdfW *= lightPickProb;
        directPdfW *= lightPickProb;

        oLightState.mThroughput /= emissionPdfW;
        oLightState.mPathLength = 1;
        oLightState.mIsFiniteLight = light->IsFinite() ? 1 : 0;

        // Light sub-path MIS quantities. Implements [tech. rep. (31)-(33)] partially.
        // The evaluation is completed after tracing the emission ray in the light sub-path loop.
        // Delta lights are handled as well [tech. rep. (48)-(50)].
        {
            oLightState.dVCM = Mis(directPdfW / emissionPdfW);

            if (!light->IsDelta())
            {
                const float usedCosLight = light->IsFinite() ? cosLight : 1.f;
                oLightState.dVC = Mis(usedCosLight / emissionPdfW);
            }
            else
            {
                oLightState.dVC = 0.f;
            }

            oLightState.dVM = oLightState.dVC * mMisVcWeightFactor;
        }
    }

    // Computes contribution of light sample to camera by splatting it onto the
    // framebuffer. Multiplies by throughput (obviously, as nothing is returned).
    void ConnectToCamera(
        const SubPathState &aLightState,
        const vec3        &aHitpoint,
        const LightBSDF    &aBsdf)
    {
        const Camera &camera = *m_scene.m_camera;
        vec3 directionToCamera = camera.position() - aHitpoint;

        // Check point is in front of camera
        if (dot(camera.forwardDirection(), -directionToCamera) <= 0.f)
        {
            return;
        }

        // Check it projects to the screen (and where)
        vec2 imagePos;
        if (!camera.rasterizePoint(aHitpoint, imagePos.x, imagePos.y))
        {
            return;
        }

        // Compute distance and normalize direction to camera
        const float distEye2 = lengthSquared(directionToCamera);
        const float distance = std::sqrt(distEye2);
        directionToCamera /= distance;

        // Get the BSDF
        float cosToCamera, bsdfDirPdfW, bsdfRevPdfW;
        const vec3 bsdfFactor = aBsdf.Evaluate(
            directionToCamera, cosToCamera, &bsdfDirPdfW, &bsdfRevPdfW);

        if (isZero(bsdfFactor))
        {
            return;
        }

        bsdfRevPdfW *= aBsdf.ContinuationProb();

        // Compute pdf conversion factor from image plane area to surface area
        const float cosAtCamera = dot(camera.forwardDirection(), -directionToCamera);
        const float imagePointToCameraDist = camera.filmDistance() / cosAtCamera;
        const float imageToSolidAngleFactor = sq(imagePointToCameraDist) / cosAtCamera;
        const float imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / sq(distance);

        // We put the virtual image plane at such a distance from the camera origin
        // that the pixel area is one and thus the image plane sampling pdf is 1.
        // The area pdf of aHitpoint as sampled from the camera is then equal to
        // the conversion factor from image plane area density to surface area density
//         const float cameraPdfA = imageToSurfaceFactor;
        const float cameraPdfA = imageToSurfaceFactor * camera.pixelPdfA();

        // Partial light sub-path weight [tech. rep. (46)]. Note the division by
        // mLightPathCount, which is the number of samples this technique uses.
        // This division also appears a few lines below in the framebuffer accumulation.
        const float wLight = Mis(cameraPdfA / mLightSubPathCount) * (
            mMisVmWeightFactor + aLightState.dVCM + aLightState.dVC * Mis(bsdfRevPdfW));

        // Partial eye sub-path weight is 0 [tech. rep. (47)]

        // Full path MIS weight [tech. rep. (37)]. No MIS for traditional light tracing.
        const float misWeight = mLightTraceOnly ? 1.0f : (1.0f / (wLight + 1.0f));

        const float surfaceToImageFactor = 1.0f / imageToSurfaceFactor;

        // We divide the contribution by surfaceToImageFactor to convert the (already
        // divided (which is always 1 ??)) pdf from surface area to image plane area, w.r.t. which the
        // pixel integral is actually defined. We also divide by the number of samples
        // this technique makes, which is equal to the number of light sub-paths
//         const vec3 contrib = misWeight * aLightState.mThroughput * bsdfFactor /
//             (mLightSubPathCount * surfaceToImageFactor);
        const float pixelArea = 1.0f / (camera.pixelPdfA()); // uniform pdf, such that pdfA * area = 1
        const float imagePdfA = surfaceToImageFactor;
        const vec3 contrib = misWeight * aLightState.mThroughput * bsdfFactor /
            (mLightSubPathCount * imagePdfA * pixelArea);

        if (!isZero(contrib))
        {
            if (m_scene.Occluded(aHitpoint, directionToCamera, distance))
            {
                return;
            }

            int pos_x = static_cast<int>(imagePos.x);
            int pos_y = static_cast<int>(imagePos.y);
            m_framebuffer.accumulatePixel(pos_x, pos_y, contrib);
            if (m_pyramid)
            {
                m_pyramid->get(aLightState.mPathLength + 1, 1)->accumulatePixel(pos_x, pos_y, contrib);
            }
        }
    }

    // Samples a scattering direction camera/light sample according to BSDF.
    // Returns false for termination
    template<bool tLightSample>
    bool SampleScattering(
        const BSDF<tLightSample> &aBsdf,
        const vec3              &aHitPoint,
        SubPathState             &aoState)
    {
        // x,y for direction, z for component. No rescaling happens
        vec3 rndTriplet = vec3(randf(), randf(), randf());
        float bsdfDirPdfW, cosThetaOut;
        uint32_t sampledEvent;

        vec3 bsdfFactor = aBsdf.Sample(rndTriplet, aoState.mDirection,
            bsdfDirPdfW, cosThetaOut, &sampledEvent);

        if (isZero(bsdfFactor))
        {
            return false;
        }

        // If we sampled specular event, then the reverse probability
        // cannot be evaluated, but we know it is exactly the same as
        // forward probability, so just set it. If non-specular event happened,
        // we evaluate the pdf
        float bsdfRevPdfW = bsdfDirPdfW;
        if ((sampledEvent & LightBSDF::kSpecular) == 0)
        {
            bsdfRevPdfW = aBsdf.Pdf(aoState.mDirection, true);
        }

        // Russian roulette
        const float contProb = aBsdf.ContinuationProb();
        if (randf() > contProb)
        {
            return false;
        }

        bsdfDirPdfW *= contProb;
        bsdfRevPdfW *= contProb;

        // Sub-path MIS quantities for the next vertex. Only partial - the
        // evaluation is completed when the actual hit point is known,
        // i.e. after tracing the ray, in the sub-path loop.

        if (sampledEvent & LightBSDF::kSpecular)
        {
            // Specular scattering case [tech. rep. (53)-(55)] (partially, as noted above)
            aoState.dVCM = 0.f;
            //aoState.dVC *= Mis(cosThetaOut / bsdfDirPdfW) * Mis(bsdfRevPdfW);
            //aoState.dVM *= Mis(cosThetaOut / bsdfDirPdfW) * Mis(bsdfRevPdfW);
            assert(bsdfDirPdfW == bsdfRevPdfW);
            aoState.dVC *= Mis(cosThetaOut);
            aoState.dVM *= Mis(cosThetaOut);

            aoState.mSpecularPath &= 1;
        }
        else
        {
            // Implements [tech. rep. (34)-(36)] (partially, as noted above)
            aoState.dVC = Mis(cosThetaOut / bsdfDirPdfW) * (
                aoState.dVC * Mis(bsdfRevPdfW) +
                aoState.dVCM + mMisVmWeightFactor);

            aoState.dVM = Mis(cosThetaOut / bsdfDirPdfW) * (
                aoState.dVM * Mis(bsdfRevPdfW) +
                aoState.dVCM * mMisVcWeightFactor + 1.f);

            aoState.dVCM = Mis(1.f / bsdfDirPdfW);

            aoState.mSpecularPath &= 0;
        }

        aoState.mOrigin = aHitPoint;
        aoState.mThroughput *= bsdfFactor * (cosThetaOut / bsdfDirPdfW);

        return true;
    }
};

void DemoVCM::run()
{
    Scene_debug scene;

//     EyeLight pt(scene);
//     PathTracer pt(scene);

//     float radiusFactor = 0.003f;
//     float radiusFactor = 0.000005f;
//     float radiusFactor = 0.000003f; // this has great influence
//     float radiusFactor = 0.000002f; // this has great influence
//     float radiusFactor = 0.000001f; // this has great influence <0>
    float radiusFactor = 0.003f; // this has great influence <1>
    float radiusAlpha = 0.75f;
//     VertexCM pt(scene, VertexCM::kLightTrace, radiusFactor, radiusAlpha);
//     VertexCM pt(scene, VertexCM::kPpm, radiusFactor, radiusAlpha);
//     VertexCM pt(scene, VertexCM::kBpm, radiusFactor, radiusAlpha);
    VertexCM pt(scene, VertexCM::kBpt, radiusFactor, radiusAlpha);
//     VertexCM pt(scene, VertexCM::kVcm, radiusFactor, radiusAlpha);

    for (int i = 0; i < 100; i++)
    {
        std::cout << i << std::endl;
        pt.RunIteration(i);
    }

    pt.dumpPyramid();

    FrameBuffer tmp;
    pt.GetFramebuffer(tmp);

    tmp.dumpHDR("www.hdr");
    tmp.tonemapGamma(2.2f);
    tmp.dumpPPM("www.ppm");
}

int runTest(int argc, char *argv[])
{
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());
    omp_set_num_threads(std::max(1, omp_get_max_threads() - 1));
    //     omp_set_num_threads(1); // single thread
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());


    DemoVCM *demo = new DemoVCM;

    demo->run();
    delete demo;

    return 0;
}

int main(int argc, char **argv)
{
    runTest(argc, argv);
}

