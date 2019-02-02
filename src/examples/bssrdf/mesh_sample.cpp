#include <cassert>
#include <algorithm>
#include "mesh_sample.h"
#include "constants.h"
#include "numeric.h"

#define DEBUGGING 0
#define LOGGING 1

typedef __int64 Integer;

struct comparatorPhaseGroupID {
    bool operator() (const std::pair<Integer, int>& a, const std::pair<Integer, int>& b) {
        return a.second < b.second;
    }
};

template <typename T1, typename T2>
struct comparatorPairSecond {
    bool operator() (const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) {
        return a.second < b.second;
    }
};

template <typename T1, typename T2>
struct comparatorPairFirst {
    bool operator() (const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) {
        return a.first < b.first;
    }
};

struct Sample {
    float x, y, z;
    Integer id; // cell id
    int pgid; // phase group id
    Integer index; // index in the input array, in order to access payloads, etc.

    struct comparator {
        bool operator() (const Sample& a, const Sample& b) const {
            return a.id < b.id;
        }
    };
    struct comparatorPhaseGroup {
        bool operator() (const Sample& a, const Sample& b) const {
            return a.pgid < b.pgid;
        }
    };
    struct printer {
        void operator() (const Sample& a) const {
            printf("[%lld / %d] %f, %f, %f\n", a.id, a.pgid, a.x, a.y, a.z);
        }
    };
    static float distanceSquared(const Sample& a, const Sample& b) {
        float x = a.x - b.x;
        float y = a.y - b.y;
        float z = a.z - b.z;
        return x*x + y*y + z*z;
    }
};

struct HashGrid {
    float m_min[3];
    float m_max[3];
    float m_delta; //cell size

    Integer m_dim[3];
    int m_perm[27]; // phase group permutation

    HashGrid(float radius) {
        m_delta = radius / std::sqrtf(3.0f);
        clear();

        int temp[27] = { 5, 6, 17, 22, 0, 23, 13, 3, 19, 2, 16, 25, 9, 24, 15, 8, 20, 7, 14, 11, 10, 21, 4, 18, 12, 1, 26 };
        for (int n = 0; n < 27; n++) {
            m_perm[n] = temp[n];
        }

#if LOGGING
        printf("delta: %f\n", m_delta);
#endif
    }
    void clear() {
        m_min[0] = m_min[1] = m_min[2] = NUM_INFINITY;
        m_max[0] = m_max[1] = m_max[2] = -NUM_INFINITY;
    }
    void computeBoundingBox(std::vector<Sample>& samples) {
        for (Integer n = 0; n < samples.size(); n++) {
            m_min[0] = f_min(m_min[0], samples[n].x);
            m_min[1] = f_min(m_min[1], samples[n].y);
            m_min[2] = f_min(m_min[2], samples[n].z);

            m_max[0] = f_max(m_max[0], samples[n].x);
            m_max[1] = f_max(m_max[1], samples[n].y);
            m_max[2] = f_max(m_max[2], samples[n].z);
        }

        for (int n = 0; n < 3; n++)
        {
            m_dim[n] = Integer(ceil((m_max[n] - m_min[n]) / m_delta));
        }

#if LOGGING
        printf("bounding box:\n");
        printf("min: %f, %f, %f\n", m_min[0], m_min[1], m_min[2]);
        printf("max: %f, %f, %f\n", m_max[0], m_max[1], m_max[2]);
        printf("dimension:\n");
        printf("%lld x %lld x %lld\n", m_dim[0], m_dim[1], m_dim[2]);
#endif

        //unit test
        //         for (int test = 0; test < 1000000; test++) {
        //             Integer ijk[3] = { rand() % m_dim[0], rand() % m_dim[1], rand() % m_dim[2] };
        //             Integer temp = coordToIndex(ijk[0], ijk[1], ijk[2]);
        //             Integer IJK[3];
        //             indexToCoord(temp, IJK);
        //             assert((ijk[0] == IJK[0]) && (ijk[1] == IJK[1]) && (ijk[2] == IJK[2]));
        //         }
    }

    Integer coordToIndex(Integer i, Integer j, Integer k) {
        return i + j * m_dim[0] + k * m_dim[0] * m_dim[1];
    }
    void indexToCoord(Integer index, Integer coords[3]) {
        Integer N = m_dim[0];
        Integer NN = m_dim[0] * m_dim[1];
        coords[2] = index / NN;
        coords[1] = (index % NN) / N;
        coords[0] = (index % NN) % N;
    }

    Integer computeCellID(const Sample& s) {
        Integer i = (s.x - m_min[0]) / m_delta;
        Integer j = (s.y - m_min[1]) / m_delta;
        Integer k = (s.z - m_min[2]) / m_delta;

        return coordToIndex(i, j, k);
    }
    int computePhaseGroupID(const Sample& s) {
        Integer i = (s.x - m_min[0]) / m_delta;
        Integer j = (s.y - m_min[1]) / m_delta;
        Integer k = (s.z - m_min[2]) / m_delta;

        return m_perm[(i % 3) + (j % 3) * 3 + (k % 3) * 9];
    }

    bool isValid(Integer i, Integer j, Integer k) {
        return (
            i >= 0 && i < m_dim[0] &&
            j >= 0 && j < m_dim[1] &&
            k >= 0 && k < m_dim[2]);

    }
};

struct HashTable {
    struct HashBucket {
        struct HashSlot {
            Integer key; // global cell ID
            Integer first_index; // indexing the samples array
            Integer sample; // indexing the samples array
        };

        std::vector< HashSlot > m_slots;
        int total;

        HashBucket() {
            m_slots.clear();
            total = 0;
        }
        void push(const HashSlot& slot) {
            m_slots.push_back(slot);
            total = m_slots.size();
        }

        struct printer {
            void operator() (const HashBucket& a) const {
                if (a.total)
                {
                    printf("[total %d] ", a.total);
                    for (int n = 0; n < a.total; n++) {
                        printf("{ %lld, %lld, %lld }, ",
                            a.m_slots[n].key, a.m_slots[n].first_index, a.m_slots[n].sample);
                    }
                    printf("\n");
                }
            }
        };
    };

    typedef HashBucket::HashSlot CellType;

    std::vector<HashBucket> m_hashTable;
    Integer m_hashSize;

    Integer hashFunc(Integer key) { return key % m_hashSize; }
    CellType *findCell(Integer cell_id) {
        Integer hash_idx = hashFunc(cell_id);
        HashBucket& bucket = m_hashTable[hash_idx];

        CellType *cell = NULL;
        for (int m = 0; m < bucket.total; m++) {
            if (bucket.m_slots[m].key == cell_id) {
                cell = &bucket.m_slots[m];
                break;
            }
        }
        return cell;
    }
    void push(const CellType& temp) {
        Integer hash_idx = hashFunc(temp.key);
        m_hashTable[hash_idx].push(temp);
    }
    void print() {
        std::for_each(m_hashTable.begin(), m_hashTable.end(), HashBucket::printer());
    }

    HashTable(Integer N) {
        m_hashTable.resize(N);
        m_hashSize = m_hashTable.size();
    }

    Integer size() { return m_hashSize; }
};

struct PointCloud {
    float m_radius;
    float m_radiusSquared;
    std::vector<Sample> m_array;
    std::vector<Position> m_output;
    Integer m_pcount;

    HashGrid *m_hashGrid;
    HashTable *m_hashTable;

    void printHashTable() {
        m_hashTable->print();
    }

    std::vector<std::pair<Integer, int>> m_validCells; // mapping between cell id and phase group id
    int phaseGroupBegin[27];
    int phaseGroupEnd[27];

    const std::vector<Position> *m_refInput;

    PointCloud(const std::vector<Position>& inputArray, float radius) {
        Integer N = inputArray.size();
        m_refInput = &inputArray;

#if LOGGING
        printf("\n\nPoisson Sampling Info:\n");
        printf("input count: %lld\nradius: %f\n", N, radius);
#endif

        m_radius = radius;
        m_radiusSquared = m_radius * m_radius;
        m_array.resize(N);
        m_pcount = m_array.size();
        m_hashGrid = new HashGrid(m_radius);
    }
    ~PointCloud() {
        delete m_hashGrid;
        delete m_hashTable;
    }

    size_t size() { return m_array.size(); }
    Sample& operator[] (size_t n) { return m_array[n]; }
    std::vector<Sample>::iterator begin() { return m_array.begin(); }
    std::vector<Sample>::iterator end() { return m_array.end(); }

    void print() {
        std::for_each(m_array.begin(), m_array.end(), Sample::printer());
    }
    void buildHashGrid() {
        m_hashGrid->computeBoundingBox(m_array);

        // sort along cell id
        for (Integer n = 0; n < m_array.size(); n++) {
            m_array[n].id = m_hashGrid->computeCellID(m_array[n]);
        }
        std::sort(m_array.begin(), m_array.end(), Sample::comparator()); // the samples array (m_array) should be sorted only once

        {
            int totalCells = 0;
            Integer lastCellID = -1;

            for (Integer n = 0; n < m_array.size(); n++) {
                if (m_array[n].id != lastCellID) {
                    totalCells++;
                }
                lastCellID = m_array[n].id;
            }

            m_hashTable = new HashTable(totalCells * 4);

#if LOGGING
            printf("total occupied cells: %lld\n", totalCells);
            printf("hash table size: %d\n", m_hashTable->size());
#endif

            m_validCells.resize(totalCells);
        }

        {
            Integer lastCellID = -1;

            for (Integer n = 0; n < m_array.size(); n++) {

                // if at the start of a cell
                if (m_array[n].id != lastCellID) {

                    HashTable::CellType temp;
                    temp.key = m_array[n].id;   // cell id
                    temp.first_index = n;   // array offset
                    temp.sample = -1;

                    m_hashTable->push(temp);
                }
                lastCellID = m_array[n].id;
            }
        }
    }

    void blueSample() {
        for (int n = 0; n < m_array.size(); n++) {
            m_array[n].pgid = m_hashGrid->computePhaseGroupID(m_array[n]);
        }

        // no more sorting on the samples array
        //         std::sort(m_array.begin(), m_array.end(), Sample::comparatorPhaseGroup());

        {
            Integer lastCellID = -1;

            int count = 0;

            // scan again over all the samples, certainly not the fastest approach here
            for (int n = 0; n < m_array.size(); n++) {

                Integer thisCellID = m_array[n].id;

                // if start of a cell
                if (thisCellID != lastCellID) {

                    m_validCells[count++] = std::make_pair(thisCellID, m_array[n].pgid);
                }
                lastCellID = thisCellID;
            }

#if DEBUGGING
            for (int n = 0; n < m_validCells.size(); n++) {
                printf("cell id [%lld], phase group id [%d]\n", m_validCells[n].first, m_validCells[n].second);
            }
            printf("=====================================================\n");
            system("pause");
#endif

            //             std::sort(m_validCells.begin(), m_validCells.end(), comparatorPhaseGroupID());
            std::sort(m_validCells.begin(), m_validCells.end(), comparatorPairSecond<Integer, int>()); // sort against phase group id

#if DEBUGGING
            for (int n = 0; n < m_validCells.size(); n++) {
                printf("cell id [%lld], phase group id [%d]\n", m_validCells[n].first, m_validCells[n].second);
            }
            printf("=====================================================\n");
            system("pause");
#endif
        }

        {
#if 0 // bug
            phaseGroupBegin[0] = 0;
            for (int n = 1; n < m_validCells.size(); n++) {
                if (m_validCells[n].second != m_validCells[n - 1].second) {
                    phaseGroupBegin[m_validCells[n].second] = n;
                    phaseGroupEnd[m_validCells[n - 1].second] = n;
                }
            }
            phaseGroupEnd[26] = m_validCells.size();
#else
            // fixed bug of empty phase group(s) in case of too few candidate points
            for (int n = 0; n < 27; n++) {
                phaseGroupBegin[n] = -1;
                phaseGroupEnd[n] = -1;
            }

#if 1 //this is (NOT) stupid
            // the first valid cell (hence the first valid phase group) may not be the first phase group
            phaseGroupBegin[m_validCells[0].second] = 0;

            for (int n = 1; n < m_validCells.size(); n++) {
                if (m_validCells[n].second != m_validCells[n - 1].second) {
                    phaseGroupBegin[m_validCells[n].second] = n;
                    phaseGroupEnd[m_validCells[n - 1].second] = n;
                }
            }

            // the last valid cell (hence the last valid phase group) may not be the last phase group
            phaseGroupEnd[m_validCells[m_validCells.size() - 1].second] = m_validCells.size();
#else // this causes seg fault
            // the first valid cell (hence the first valid phase group) may not be the first phase group
            // the last valid cell (hence the last valid phase group) may not be the last phase group
            int total = m_validCells.size();
            for (int n = 0; n <= total; n++) {
                if (m_validCells[n].second != m_validCells[n - 1].second) {
                    if (n < total) phaseGroupBegin[m_validCells[n].second] = n;
                    if (n > 0) phaseGroupEnd[m_validCells[n - 1].second] = n;
                }
            }
#endif
#endif
        }

#if DEBUGGING
        for (int n = 0; n < 27; n++) {
            printf("phase group %d = cell range [%d, %d]\n", n, phaseGroupBegin[n], phaseGroupEnd[n]);
        }
        system("pause");
#endif

        for (int t = 0; t < 30; t++)
        {
            //for each phase group
            for (int n = 0; n < 27; n++) {

                int begin = phaseGroupBegin[n];
                int end = phaseGroupEnd[n];

                // parallel for each cell in phase group
#pragma omp parallel for
                for (int i = begin; i < end; i++) {

                    Integer cell_id = m_validCells[i].first;
                    HashTable::CellType *cell = m_hashTable->findCell(cell_id);

                    //                     if (cell->sample >= 0) {
                    //                         //                         printf("occupied by sample %d\n", cell->sample);
                    //                         continue;
                    //                     }

                    // debugging for Release
                    // if (cell->key != cell_id) { 
                    //     printf("cell->key(%d) != cell_id(%d)\n", cell->key, cell_id);
                    //     while(1);
                    // }

                    assert(cell->key == cell_id);

#if DEBUGGING
                    printf("cell = %lld, %lld, %lld\n", cell->key, cell->first_index, cell->sample);
#endif

                    Integer randomSampleIndex = cell->first_index + t;
                    randomSampleIndex = randomSampleIndex % m_pcount; // trick to avoid range overflow
                    Sample& randomSample = m_array[randomSampleIndex];

                    if (cell->sample == randomSampleIndex) {
                        continue;
                    }

                    if (randomSample.id == cell_id) // random sample in this cell
                    {
                        // debugging for Release
                        // if (randomSample.pgid != n) { 
                        //     printf("randomSample.pgid(%d) != n(%d)\n", randomSample.pgid, n);
                        //     printf("\tbegin(%d), end(%d)\n", begin, end);
                        //     while(1);
                        // } // there must have been too few candidates

                        assert(randomSample.pgid == n);

                        Integer ijk[3];
                        m_hashGrid->indexToCoord(cell_id, ijk);

#if DEBUGGING
                        printf("ijk = %lld, %lld, %lld\n", ijk[0], ijk[1], ijk[2]);
                        printf("neighbors:\n\t");
#endif
                        bool conflict = false;
                        for (int K = -2; K <= 2; K++) {            // FIXED BUG : must search 2-neighbor not only 1-neighbor
                            for (int J = -2; J <= 2; J++) {        // FIXED BUG : must search 2-neighbor not only 1-neighbor
                                for (int I = -2; I <= 2; I++) {    // FIXED BUG : must search 2-neighbor not only 1-neighbor
                                    Integer ii = ijk[0] + I;
                                    Integer jj = ijk[1] + J;
                                    Integer kk = ijk[2] + K; // J; !!!!!!!!! FIXED TYPO !!!!!!!!!!!
                                    if ((I != 0 || J != 0 || K != 0) && m_hashGrid->isValid(ii, jj, kk)) // FIXED BUG : do not check self !!!!!
                                    {
                                        Integer neighbor_cell_id = m_hashGrid->coordToIndex(ii, jj, kk);
#if DEBUGGING
                                        printf("nid(%lld), ", neighbor_cell_id);
#endif

                                        HashTable::CellType *neighbor_cell = m_hashTable->findCell(neighbor_cell_id);
                                        if (neighbor_cell && neighbor_cell->sample >= 0) { // it is important because the distribution may be sparse, and 
                                            // neighbor cells don't exist at all

#if DEBUGGING
                                            printf("nref(%lld), ", neighbor_cell->sample);
#endif

                                            if (Sample::distanceSquared(randomSample, m_array[neighbor_cell->sample]) < m_radiusSquared) {
                                                conflict = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }

#if DEBUGGING
                        printf("\n");
#endif

                        if (!conflict) {
                            cell->sample = randomSampleIndex;
                        }
                    }
                }
            }
        }

#if 0
        for (Integer n = 0; n < m_hashTable.size(); n++) {
            HashBucket& refBucket = m_hashTable[n];
            for (int m = 0; m < refBucket.total; m++) {
                HashBucket::HashSlot& refSlot = refBucket.m_slots[m];
                Integer sample_id = refSlot.sample;
                if (sample_id >= 0) {
                    printf("cell %lld : sample %lld\n", refSlot.key, sample_id);
                }
            }
        }
#else
        std::sort(m_validCells.begin(), m_validCells.end(), comparatorPairFirst<Integer, int>()); // sort against cell id

        for (Integer n = 0; n < m_validCells.size(); n++) {
            Integer cell_id = m_validCells[n].first;
            HashTable::CellType *cell = m_hashTable->findCell(cell_id);
            Integer sample = cell->sample;
            if (sample >= 0) {
#if DEBUGGING
                printf("cell %lld : sample %lld\n", cell->key, sample);
#endif
                //                 m_output.push_back(Position(m_array[sample].x, m_array[sample].y, m_array[sample].z)); // not considering payload!!!!!
                //                 m_output.push_back(Position(m_array[sample].x, m_array[sample].y, m_array[sample].z, (*m_refInput)[sample].payload));
                //                 m_output.push_back((*m_refInput)[sample]); // does need to not involve payloads, etc. // not correct
                m_output.push_back((*m_refInput)[m_array[sample].index]); // does need to not involve payloads, etc., index is for retro-indexing
            }
        }
#endif

#if LOGGING
        printf("output count: %lld\n", m_output.size());
#endif
    }
};

int poissonSample(const std::vector<Position>& p_in, std::vector<Position>& p_out, float min_dist) {

    //     PointCloud p(1000, 0.10005);
    //     for (Integer n = 0; n < p.size(); n++) {
    //         Sample temp;
    //         temp.x = randf();
    //         temp.y = randf();
    //         temp.z = randf();
    //         temp.id = -1;
    //         temp.pgid = -1;
    // 
    //         p[n] = temp;
    //     }

    PointCloud p(p_in, min_dist);

    for (Integer n = 0; n < p.size(); n++) {
        Sample temp;
        temp.x = p_in[n].x;
        temp.y = p_in[n].y;
        temp.z = p_in[n].z;
        temp.id = -1;
        temp.pgid = -1;
        temp.index = n;

        p[n] = temp;
    }

    p.buildHashGrid();

#if DEBUGGING
    printf("=====================================================\n");
    system("pause");

    p.print();
    printf("=====================================================\n");
    system("pause");

    p.printHashTable();
    printf("=====================================================\n");
    system("pause");
#endif

    p.blueSample();

#if DEBUGGING
    printf("=====================================================\n");
    system("pause");

    p.print();
    printf("=====================================================\n");
    system("pause");

    p.printHashTable();
    printf("=====================================================\n");
    system("pause");
#endif

    //     printf("=====================================================\n");
    //     system("pause");

    p_out = p.m_output;

    return 0;
}



void check_min_distance(std::vector<Position>& pp)
{
    float mindist = NUM_INFINITY;
    for (int i = 0; i < pp.size() - 1; i++) {
        for (int j = i + 1; j < pp.size(); j++) {
            float dist = Position::distanceSquared(pp[i], pp[j]);
            mindist = f_min(mindist, dist);
        }
    }
    printf("mindist of output: %f\n", std::sqrtf(mindist));
}



//         m_hashTable[0].push(m_array[0].id);
//         for (Integer n = 1; n < m_array.size(); n++) {
//             if (m_array[n].id != m_array[n-1].id) {
//                 Integer hash_idx = m_array[n].id % m_hashTable.size();
//                 m_hashTable[hash_idx].push(m_array[n].id);
//             }
//         }


float Position::distanceSquared(const Position& a, const Position& b)
{
    float x = a.x - b.x;
    float y = a.y - b.y;
    float z = a.z - b.z;
    return x*x + y*y + z*z;
}

Position::Position(float x, float y, float z) : x(x), y(y), z(z)
{

}

Position::Position()
{

}
