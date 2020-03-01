//
// Created by cothrax on 3/1/20.
//

#include "stat.h"
#include "filter.h"
#include<fstream>
#include<iostream>
#include <cstring>

const char *match_opt_fn = "matches.txt";
const char *hit_opt_fn = "hits.txt";
const char *match_stat_fn = "matches_stat.dat";
const char *his_stat_fn = "hits_stat.dat";
//const char *fig_fn = "stat.jpg";

void convert_all()
{
    DataHolder ref, reads;
    reads.load("/home/cothrax/data/MH0006.5000.fa");
    cout << "reads:\t" << reads.seq_num << " " << reads.ofs[reads.seq_num] << endl;

    ref.load("/home/cothrax/data/humangut.194.ref.genomes.fa");
    cout << "ref:\t" << ref.seq_num << " " << ref.ofs[ref.seq_num] << endl;

    StatConverter match_converter("/home/cothrax/data/matches.dat");
    StatConverter hit_converter("/home/cothrax/data/hits.dat");
//
//    bit64_t tmp1, tmp2;
//    match_converter.next_pair(tmp1, tmp2);


    match_converter.convert_hollow(ref, reads, match_opt_fn);
    hit_converter.convert_hollow(ref, reads, hit_opt_fn);
//    match_converter.convert(ref, reads, match_opt_fn);
//    hit_converter.convert(ref, reads, hit_opt_fn);
}

void _plot_stat(LSHFilter &filter, const char *opt_fn, const char *stat_fn)
{
    using namespace std;

//    LSHFilter filter;
//    filter.load_projection("/media/cothrax/Elements/proj.dat");

    ifstream fin(opt_fn);
    ofstream fout(stat_fn, ios_base::out | ios_base::binary);
    char s1[100], s2[100];
    while(fin >> s1 >> s2)
    {
        lshv_t val1 = filter.lsh(s1, strlen(s1));
        lshv_t val2 = filter.lsh(s2, strlen(s2));
//        fout << __builtin_popcount(val1 ^ val2) << " ";
        char cnt = __builtin_popcountll(val1 ^ val2);
//        printf("%s %s %u %u\n", s1, s2, val1, val2);
        fout.write(&cnt, sizeof(char));
    }
    fin.close();
    fout.close();
}

void test_filter(LSHFilter &filter)
{
    _plot_stat(filter, match_opt_fn, match_stat_fn);
    _plot_stat(filter, hit_opt_fn, his_stat_fn);
    char cmd[100];
    sprintf(cmd, "python plot_stat.py %s %s", match_stat_fn, his_stat_fn);
    cout << cmd << '\n';
    system(cmd);
}

void _plot_stat_hollow(LSHFilter &filter, const char *opt_fn, const char *stat_fn)
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
        lshv_t val1 = filter.hollow_lsh(s1, strlen(s1), begin1);
        lshv_t val2 = filter.hollow_lsh(s2, strlen(s2), begin2);
//        fout << __builtin_popcount(val1 ^ val2) << " ";
        char cnt = __builtin_popcountll(val1 ^ val2);
//        printf("%s %s %u %u\n", s1, s2, val1, val2);
        fout.write(&cnt, sizeof(char));
    }
    fin.close();
    fout.close();
}

void test_filter_hollow(LSHFilter &filter)
{
    _plot_stat_hollow(filter, match_opt_fn, match_stat_fn);
    _plot_stat_hollow(filter, hit_opt_fn, his_stat_fn);
    char cmd[100];
    sprintf(cmd, "python plot_stat.py %s %s", match_stat_fn, his_stat_fn);
    cout << cmd << '\n';
    system(cmd);
}

int main()
{
    convert_all();
//    plot_stat();
    LSHFilter filter;
    RandomProjectionGenerator generator(2, 32, 6);
    filter.load_projection(generator);
//    filter.load_projection("/media/cothrax/Elements/proj.dat");
    test_filter_hollow(filter);
}