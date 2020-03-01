//
// Created by cothrax on 2/27/20.
//

#ifndef LIBJPEG_FILTER_H
#define LIBJPEG_FILTER_H

#include "param.h"
typedef bit64_t lshv_t;

class ProjectionGenerator
{
protected:
    void compress_vector();

public:
    int _k, _nbits, total;
    int **vec, **pos, **neg;

    ProjectionGenerator(int k, int nbits);
    virtual ~ProjectionGenerator();
    virtual void generate() = 0;
    virtual void dump(const char *filename);
};

class RandomProjectionGenerator: public ProjectionGenerator
{
public:
    int _density;
    RandomProjectionGenerator(int k, int nbits, int density): ProjectionGenerator(k, nbits), _density(density) {}
    void rand_proj(int seed, int *a_vec);
    virtual void  generate();
};

class LSHFilter
{
protected:
    lshv_t dict[128];
    int **pos, **neg;
    inline void update(int *acc, int key, lshv_t &val)
    {
        for(int j = 1; j <= pos[key][0]; j++)
            if(++acc[pos[key][j]] == 0) val ^= (lshv_t)1<<pos[key][j];
        for(int j = 1; j <= neg[key][0]; j++)
            if(acc[neg[key][j]]-- == 0) val ^= (lshv_t)1<<neg[key][j];
    }
    inline void update_inverse(int *acc, int key, lshv_t &val)
    {
        for(int j = 1; j <= pos[key][0]; j++)
            if(acc[pos[key][j]]-- == 0) val ^= (lshv_t)1<<pos[key][j];
        for(int j = 1; j <= neg[key][0]; j++)
            if(++acc[neg[key][j]] == 0) val ^= (lshv_t)1<<neg[key][j];
    }

public:
    int k, nbits, total, nkmer, _seed_size, frag_len;
    explicit LSHFilter(int seed_size = 11);
    ~LSHFilter();
    void load_projection(const char *filename);
    void load_projection(ProjectionGenerator &generator);
    lshv_t lsh(char *seq, int seq_len);
    lshv_t hollow_lsh(char *seq, int seq_len, int begin);
    void calc_all_lsh(char *seq, int seq_len, lshv_t *res);
    void calc_seq_lsh(char *seq, int seq_len, int seed_step, lshv_t *res);
};

class RefLSHHolder
{

};

#endif //LIBJPEG_FILTER_H
