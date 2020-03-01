//
// Created by cothrax on 2/27/20.
//

#include "filter.h"
#include<fstream>
#include<cstring>
#include<iostream>
#include<cstdio>
#include<random>
#include<cmath>

// -------------------- LSHFilter --------------------

LSHFilter::LSHFilter(int seed_size)
{
    _seed_size = seed_size;
    frag_len = seed_size * 3;
    pos = nullptr;
    neg = nullptr;

    memset(dict, 0, sizeof(dict));
    dict[1] = dict['C'] = dict['c'] = 1;
    dict[2] = dict['G'] = dict['g'] = 2;
    dict[3] = dict['T'] = dict['t'] = 3;
}

LSHFilter::~LSHFilter()
{
    if(pos)
    {
        for(int i = 0; i < total; i++) delete [] pos[i];
        delete [] pos;
    }
    if(neg)
    {
        for(int i = 0; i < total; i++) delete [] neg[i];
        delete [] neg;
    }
}

void LSHFilter::load_projection(const char *filename)
{
    using std::ifstream;
    using std::ios_base;
    ifstream fin(filename, ios_base::in | ios_base::binary);
    fin.read((char *) &nbits, sizeof(nbits));
    fin.read((char *) &k, sizeof(k));
    total = 1<<(k<<1);
    nkmer = frag_len - k + 1;

    pos = new int * [total];
    neg = new int * [total];
    for(int i = 0; i < total; i++)
    {
        int tmp;
        fin.read((char *) &tmp, sizeof(tmp));
        pos[i] = new int[tmp+1];
        pos[i][0] = tmp;
        for(int j = 1; j <= tmp; j++) fin.read((char *) &pos[i][j], sizeof(pos[i][j]));

        fin.read((char *) &tmp, sizeof(tmp));
        neg[i] = new int[tmp+1];
        neg[i][0] = tmp;
        for(int j = 1; j <= tmp; j++) fin.read((char *) &neg[i][j], sizeof(neg[i][j]));
    }
    fin.close();
}

void LSHFilter::load_projection(ProjectionGenerator &generator)
{
    char filename[100];
    sprintf(filename, "tmp%d.dat", int(time(0)));
    generator.dump(filename);
    load_projection(filename);
    remove(filename);
}

lshv_t LSHFilter::lsh(char *seq, int seq_len)
{
    int len = strlen(seq);
    int key = 0;
    for(int j = 0; j < k; j++) key = (key<<2) | dict[seq[j]];

    int acc[nbits];
    memset(acc, 0, sizeof(acc));
    for(int i = k-1; i < len; i++)
    {
        for(int j = 1; j <= pos[key][0]; j++) acc[pos[key][j]]++;
        for(int j = 1; j <= neg[key][0]; j++) acc[neg[key][j]]--;
        key -= dict[seq[i-k+1]]<<((k-1)<<1);
        key = (key<<2) | dict[seq[i+1]];
    }

    lshv_t ret = 0;
    for(int i = 0; i < nbits; i++) ret |= (lshv_t)(acc[i] >= 0) << i;
    return ret;
}


void LSHFilter::calc_all_lsh(char *seq, int seq_len, lshv_t *res)
{
    int key = 0, rbound = seq_len - (_seed_size<<1), p;
    int vec[nkmer], acc[nbits];

    memset(acc, 0, sizeof(acc));
    for(int j = 0; j < k; j++) key = (key<<2) | dict[seq[j]];
    // memset(vec, -1, sizeof(vec));
    // printf("nkmer %d\n", nkmer);

    res[k-2] = ((bit64_t)1<<nbits)-1;
    for(int i = k-1; i < seq_len; i++)
    {
        p = i%nkmer;
        res[i] = res[i-1];
        if(i >= frag_len) update_inverse(acc, vec[p], res[i]);
        update(acc, vec[p] = key, res[i]);
        key -= dict[seq[i-k+1]]<<((k-1)<<1);
        key = (key<<2) | dict[seq[i+1]];
    }


    for(int i = 0; i < seq_len; i++)
        if(i < _seed_size) res[i] = res[frag_len-1];
        else if(i > rbound) res[i] = res[rbound];
        else res[i] = res[i+(_seed_size<<1)-1];
}

void LSHFilter::calc_seq_lsh(char *seq, int seq_len, int seed_step, lshv_t *res)
{
    int key = 0, cnt = 0, rbound = seq_len - (_seed_size<<1);
    lshv_t val = ((bit64_t)1<<nbits)-1;
    int vec[nkmer], acc[nbits];
    int p, end;

    memset(acc, 0, sizeof(acc));
    for(int j = 0; j < k; j++) key = (key<<2) | dict[seq[j]];
    for(int i = k-1; i < frag_len; i++)
    {
        p = i%nkmer;
        update(acc, vec[p] = key, val);
        key -= dict[seq[i-k+1]]<<((k-1)<<1);
        key = (key<<2) | dict[seq[i+1]];
    }

    res[cnt++] = val;
    for(int i = 1; i < seq_len; i++)
    {
        if(i > _seed_size && i <= rbound)
        {
            end = i + (_seed_size << 1) - 1;
            p = end%nkmer;

            update_inverse(acc, vec[p], val);
            update(acc, vec[p] = key, val);

            // i = _seed_size+1 -> end = frag_len
            // end = i - _seed_size - 1 + frag_len = i + (_seed_size << 1) - 1
            key -= dict[seq[end-k+1]]<<((k-1)<<1);
            key = (key<<2) | dict[seq[end+1]];
        }
        if(i % seed_step == 0) res[cnt++] = val;
    }

    if(seq_len % seed_step) res[cnt++] = val;
}

lshv_t LSHFilter::hollow_lsh(char *seq, int seq_len, int begin)
{
    using std::max;
    using std::min;
    int acc[nbits];
    memset(acc, 0, sizeof(acc));

    int range[2][2] =
            {
                    {0, max(k, begin)},
                    {min(seq_len-k, begin+_seed_size), seq_len}
            };

    for(int q = 0; q < 2; q++)
    {
        int key = 0;
        for(int j = range[q][0]; j < range[q][0] + k; j++) key = (key<<2) | dict[seq[j]];
        for(int i = range[q][0] + k - 1; i < range[q][1]; i++)
        {
            for(int j = 1; j <= pos[key][0]; j++) acc[pos[key][j]]++;
            for(int j = 1; j <= neg[key][0]; j++) acc[neg[key][j]]--;
            key -= dict[seq[i-k+1]]<<((k-1)<<1);
            key = (key<<2) | dict[seq[i+1]];
        }
    }

    lshv_t ret = 0;
    for(int i = 0; i < nbits; i++) ret |= (lshv_t)(acc[i] >= 0) << i;
    return ret;
}

// -------------------- ProjectionGenerator --------------------
ProjectionGenerator::ProjectionGenerator(int k, int nbits)
{
    _k = k;
    _nbits = nbits;
    total = 1<<(k<<1);
    vec = pos = neg = nullptr;
}

ProjectionGenerator::~ProjectionGenerator()
{
    int ** tmp[3] = {vec, pos, neg};
    for(int j = 0; j < 3; j++)
        if(tmp[j])
        {
            for(int i = 0; i < total; i++) delete [] tmp[j][i];
            delete [] tmp[j];
        }
}

void ProjectionGenerator::dump(const char *filename)
{
    if(!pos || !neg) compress_vector();

    using std::ofstream;
    using std::ios_base;
    ofstream fout(filename, ios_base::binary | ios_base::out);
    fout.write((char *) &_nbits, sizeof(_nbits));
    fout.write((char *) &_k, sizeof(_k));
    for(int i = 0; i < total; i++)
    {
        for(int j = 0; j <= pos[i][0]; j++) fout.write((char *) &pos[i][j], sizeof(pos[i][j]));
        for(int j = 0; j <= neg[i][0]; j++) fout.write((char *) &neg[i][j], sizeof(neg[i][j]));
    }
    fout.close();
}

void ProjectionGenerator::compress_vector()
{
    if(!vec) generate();
//    if(pos) delete [] pos;
//    if(neg) delete [] neg;
    pos = new int * [total];
    neg = new int * [total];

    for(int i = 0; i < total; i++)
    {
        int count[3] = {0, 0, 0};
        for(int j = 0; j < _nbits; j++) count[vec[i][j]+1]++;
        pos[i] = new int[count[2]+1];
        neg[i] = new int[count[0]+1];
        pos[i][0] = neg[i][0] = 0;
        for(int j = 0; j < _nbits; j++)
            if(vec[i][j] == -1) neg[i][++neg[i][0]] = j;
            else if(vec[i][j] == 1) pos[i][++pos[i][0]] = j;
    }
}

// -------------------- RandomProjectionGenerator --------------------
void RandomProjectionGenerator::rand_proj(int seed, int *a_vec)
{
    using std::mt19937;
    mt19937 rng32(seed);
    for (int i = 0; i < _nbits; i++) {
        int cur = rng32() % _density;
        if (cur == 0) a_vec[i] = 1;
        else if (cur == 1) a_vec[i] = -1;
        else a_vec[i] = 0;
    }
}

void RandomProjectionGenerator::generate()
{
    vec = new int * [total];
    for(int i = 0; i < total; i++)
    {
        vec[i] = new int [_nbits];
        rand_proj(i, vec[i]);
    }
}

/*
inline void debug()
{
    using namespace std;
    LSHFilter filter;
    filter.load_projection("/media/cothrax/Elements/proj.dat");
    lshv_t res1[1000], res2[1000], res3[1000];
    char seq[1000];
    while(1)
    {
        cin >> seq;
        int len = strlen(seq);
        int seed_step = 10;
        filter.calc_all_lsh(seq, len, res1);
        filter.calc_seq_lsh(seq, len, 1, res2);
        filter.calc_seq_lsh(seq, len, seed_step, res3);

        for(int i = 0, j = 0; i < len; i++)
        {
            cout << res1[i] << "\t" << res2[i];
            if(i % seed_step == 0 || i + 1 == len) cout << '\t' << res3[j++];
            cout << endl;
        }
    }
}
 */
/*
void debug()
{
    using namespace std;
    RandomProjectionGenerator generator(4, 32, 6);
    generator.dump("/media/cothrax/Elements/proj_rand.dat");
    LSHFilter filter;
    filter.load_projection(generator);
    char s[1000];
    lshv_t res[1000];
    while(1)
    {
        scanf("%s", s);
        filter.calc_all_lsh(s, strlen(s), res);
        for(int i = 0; i < strlen(s); i++) cout << i << '\t' << res[i] << endl;
    }
}

int main()
{
    debug();
}
 */