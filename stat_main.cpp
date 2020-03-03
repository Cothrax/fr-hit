//
// Created by cothrax on 3/1/20.
//

#include "stat.h"
#include "filter.h"
#include<fstream>
#include<iostream>
#include<cstring>
#include<queue>
#include<algorithm>

typedef pair<int, int> Hit;
const char *match_opt_fn = "matches.txt";
const char *hit_opt_fn = "hits.txt";
const char *match_stat_fn = "matches_stat.dat";
const char *his_stat_fn = "hits_stat.dat";
//const char *fig_fn = "stat.jpg";
/*
lshv_t lsh_combined(LSHFilter &filter, char *seq, int seq_len, int begin)
{
    lshv_t val1 = filter.lsh(seq, begin);
    lshv_t val2 = filter.lsh(seq+begin+filter._seed_size, seq_len-begin-filter._seed_size);
    return (val1 << filter.nbits) | val2;
}

short evaluate(LSHFilter &filter, lshv_t val1, lshv_t val2)
{
    lshv_t msk = ((bit64_t)1<<filter.nbits) - 1;
    lshv_t res = val1 ^ val2;
    short s1 =  __builtin_popcountll(res & msk);
    short s2 = __builtin_popcountll(res & (msk<<filter.nbits));
    printf("%d %d %d\n", s1, s2, s1|s2 ? (2*s1*s2)/(s1+s2) : 0);
    return s1|s2 ? (2*s1*s2)/(s1+s2) : 0;
//    return s1 * s2;
}

void debug()
{
    DataHolder ref, reads;
    reads.load("/home/cothrax/data/MH0006.5000.fa");
    cout << "reads:\t" << reads.seq_num << " " << reads.ofs[reads.seq_num] << endl;

    ref.load("/home/cothrax/data/humangut.194.ref.genomes.fa");
    cout << "ref:\t" << ref.seq_num << " " << ref.ofs[ref.seq_num] << endl;

    StatConverter hit_converter("/home/cothrax/data/hits.dat");
    hit_converter.convert_hollow(ref, reads, hit_opt_fn);
}
*/

void convert_all(bool hollow)
{
    DataHolder ref, reads;
    reads.load("/home/cothrax/data/MH0006.5000.fa");
    cout << "reads:\t" << reads.seq_num << " " << reads.ofs[reads.seq_num] << endl;

    ref.load("/home/cothrax/data/humangut.194.ref.genomes.fa");
    cout << "ref:\t" << ref.seq_num << " " << ref.ofs[ref.seq_num] << endl;

    StatConverter match_converter("/home/cothrax/data/matches.dat");
    StatConverter hit_converter("/home/cothrax/data/hits.dat");

    if(hollow)
    {
        match_converter.convert_hollow(ref, reads, match_opt_fn);
        hit_converter.convert_hollow(ref, reads, hit_opt_fn);
    }
    else
    {
        match_converter.convert(ref, reads, match_opt_fn);
        hit_converter.convert(ref, reads, hit_opt_fn);
    }
}

void _plot_stat(LSHFilter &filter, const char *opt_fn, const char *stat_fn, bool hollow)
{
    using namespace std;

//    LSHFilter filter;
//    filter.load_projection("/media/cothrax/Elements/proj.dat");

    ifstream fin(opt_fn);
    ofstream fout(stat_fn, ios_base::out | ios_base::binary);
    char s1[100], s2[100];
    while(fin >> s1 >> s2)
    {
        if(s1[0] == '0') continue;

        lshv_t val1, val2;
        if(hollow)
        {
            int begin1, begin2;
            fin >> begin1 >> begin2;
            val1 = filter.hollow_lsh(s1, strlen(s1), begin1);
            val2 = filter.hollow_lsh(s2, strlen(s2), begin2);
        }
        else
        {
            val1 = filter.lsh(s1, strlen(s1));
            val2 = filter.lsh(s2, strlen(s2));
        }

//        fout << __builtin_popcount(val1 ^ val2) << " ";
        char cnt = __builtin_popcountll(val1 ^ val2);
//        printf("%s %s %u %u\n", s1, s2, val1, val2);
        fout.write(&cnt, sizeof(char));
    }
    fin.close();
    fout.close();
}

void print_loss(vector<Hit> &loss, const char *opt_fn, bool hollow)
{
    cout << "loss: " <<  loss.size() << endl;
    ifstream fin(opt_fn);
    sort(loss.begin(), loss.end());
    int ptr = 0;
    char s1[100], s2[100];
    for(int i = 0; fin >> s1 >> s2; )
    {
        if(s1[0] == '0') continue;
        while(ptr < loss.size() && loss[ptr].first < i) ptr++;
//        cout << i << " " << ptr << " " << loss[ptr].first << endl;
        int begin1, begin2;
        if(hollow) fin >> begin1 >> begin2;

        if(ptr < loss.size() && loss[ptr].first == i)
        {
            cout << loss[ptr].second << " " << s1 << " " << s2;
            if(hollow) cout << " " << begin1 << " " << begin2;
            cout << endl;
        }
        i++;
    }
    fin.close();
}

void _plot_stat_select(LSHFilter &filter, const char *opt_fn, const char *stat_fn, bool hollow)
{
    int limit = 7;
    using namespace std;

//    LSHFilter filter;
//    filter.load_projection("/media/cothrax/Elements/proj.dat");

    ifstream fin(opt_fn);
    ofstream fout(stat_fn, ios_base::out | ios_base::binary);
    char s1[100], s2[100];
    std::priority_queue<Hit, std::vector<Hit>, std::greater<Hit> > queue;
    vector<Hit> loss;
    for(int i = 0; fin >> s1 >> s2; )
    {
        if(s1[0] == '0')
        {
            Hit hit = queue.top();
            char tmp = hit.first;
            fout.write((char *) &tmp, sizeof(tmp));
            queue.pop();
            if(hit.first > limit) loss.push_back(make_pair(hit.second, hit.first));

            hit = queue.top();
            tmp = hit.first;
            fout.write((char *) &tmp, sizeof(tmp));
            if(hit.first > limit) loss.push_back(make_pair(hit.second, hit.first));

            while(!queue.empty()) queue.pop();
            continue;
        }

        lshv_t val1, val2;
        if(hollow)
        {
            int begin1, begin2;
            fin >> begin1 >> begin2;
            val1 = filter.hollow_lsh(s1, strlen(s1), begin1);
            val2 = filter.hollow_lsh(s2, strlen(s2), begin2);
        }
        else
        {
            val1 = filter.lsh(s1, strlen(s1));
            val2 = filter.lsh(s2, strlen(s2));
        }

//        fout << __builtin_popcount(val1 ^ val2) << " ";
        char cnt = __builtin_popcountll(val1 ^ val2);
//        printf("%s %s %u %u\n", s1, s2, val1, val2);
//        fout.write(&cnt, sizeof(char));
        queue.push(make_pair((int)cnt, i++));
    }
    fin.close();
    fout.close();

    print_loss(loss, opt_fn, hollow);
}

void test_filter(LSHFilter &filter, const char *fig_fn)
{
    char cmd[100];
    sprintf(cmd, "python plot_stat.py %s %s %s", match_stat_fn, his_stat_fn, fig_fn);
    cout << cmd << '\n';
    system(cmd);
}

/*
void _plot_stat_combined(LSHFilter &filter, const char *opt_fn, const char *stat_fn)
{
    using namespace std;

//    LSHFilter filter;
//    filter.load_projection("/media/cothrax/Elements/proj.dat");

    ifstream fin(opt_fn);
    ofstream fout(stat_fn, ios_base::out | ios_base::binary);
    char s1[100], s2[100];
    int begin1, begin2;

    while(fin >> begin1 >> begin2 >> s1 >> s2)
    {
        lshv_t val1 = lsh_combined(filter, s1, strlen(s1), begin1);
        lshv_t val2 = lsh_combined(filter, s2, strlen(s2), begin2);
//        fout << __builtin_popcount(val1 ^ val2) << " ";
        short cnt = evaluate(filter, val1, val2);
//        if(strlen(s1) < 30) continue;
//        if(cnt > 16)
//            printf("%s %s %d %d\n", s1, s2, begin1, begin2);
        fout.write((char *)&cnt, sizeof(short));
    }
    fin.close();
    fout.close();
}


void test_filter_combined(LSHFilter &filter)
{
    _plot_stat_combined(filter, match_opt_fn, match_stat_fn);
    _plot_stat_combined(filter, hit_opt_fn, his_stat_fn);
    char cmd[100];
    sprintf(cmd, "python plot_stat.py %s %s", match_stat_fn, his_stat_fn);
    cout << cmd << '\n';
    system(cmd);
}
*/
void test_random_filter(int k, int nbits, int density, bool hollow, bool need_convert = true)
{
    printf("args: k = %d,\tnbits = %d,\tdensity = %d,\thollow=%d\n", k, nbits, density, hollow);
    if(need_convert) convert_all(hollow);
    LSHFilter filter;
    RandomProjectionGenerator generator(k, nbits, density) ;
    filter.load_projection(generator);

    _plot_stat(filter, match_opt_fn, match_stat_fn, hollow);
    _plot_stat(filter, hit_opt_fn, his_stat_fn, hollow);

    char fn[1000];
    sprintf(fn, "k%d-n%d-d%d-h%d.jpg", k, nbits, density, (int)hollow);

    test_filter(filter, fn);
    putchar('\n');
}

void test_random_filter_select(int k, int nbits, int density, bool hollow, bool need_convert = true)
{
    printf("args: k = %d,\tnbits = %d,\tdensity = %d,\thollow=%d\n", k, nbits, density, hollow);
    if(need_convert) convert_all(hollow);
    LSHFilter filter;
    RandomProjectionGenerator generator(k, nbits, density) ;
    filter.load_projection(generator);

    _plot_stat_select(filter, match_opt_fn, match_stat_fn, hollow);
    _plot_stat(filter, hit_opt_fn, his_stat_fn, hollow);

    char fn[1000];
    sprintf(fn, "k%d-n%d-d%d-h%d.jpg", k, nbits, density, (int)hollow);

    test_filter(filter, fn);
    putchar('\n');
}


int main()
{
    int nbits = 32;
    int hollow = 1;
    test_random_filter_select(2, 16, 6, 1, 1);


    return 0;
    for(int k = 2; k <= 4; k++)
        for(int density = 2; density <= 8; density++)
            test_random_filter(k, nbits, density, hollow, k == 2 && density == 2);
}