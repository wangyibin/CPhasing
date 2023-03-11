#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include<ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>
#include <getopt.h>
#include<getopt.h>
#include<sys/sysinfo.h>

#include<algorithm>
#include<map>
#include<vector>
#include<string>
#include<sstream>

#include<cstdio>
#include <thread>
#include <iostream>
#include <fstream>
#include <chrono>
#include<functional>
#include<sys/time.h>
#include<condition_variable>
#include<mutex>
#include <unordered_map>

# define TMP_LEN 344
# define ORIENTION 4

using namespace std;


vector<string> s_set;
map<string, int> LenMap;


int get_max_pair(string input_file, string output_file){
    ifstream is(input_file, ios::in);
    ofstream ot(output_file, ios::out);
    if (!is.good() || !ot.good()){
        cerr<<"please input the right one!"<<endl;
        exit(-1);
    }

    int row = 0;
    map<float, pair<string, string>> OrientMap;
    map<pair<string, string>, float> Max_scoreMap;
    
    string chrome_1;
    string chrome_2;
    float score;

    while (!is.eof()){
        is>>chrome_1;
        is>>chrome_2;
        is>>score;
        OrientMap.insert(pair<float, pair<string, string>>({score, {chrome_1, chrome_2}}));
        row++;
        if (row == ORIENTION){
            auto iter=OrientMap.rbegin();
            Max_scoreMap.insert(pair<pair<string, string>, float>({{iter->second.first, iter->second.second}, iter->first}));
            row = 0;
            OrientMap.clear();
        }
    }

    for(auto iter = Max_scoreMap.begin(); iter != Max_scoreMap.end(); iter++){
        ot<<iter->first.first<<'\t'<<iter->first.second<<'\t'<<iter->second<<endl;
    }


    is.close();
    ot.close();
    return 0;
}




int get_distance(string input_file, string chrom_group){
    ifstream is(input_file, ios::in);
    if (!is.good()){
        cerr<<"the file is error, please input the right file! ! !"<<endl;
        exit(-1);
    }

    string suftmp;
    string chrom;
    int recount;
    int length;
    int mid_len;
    int ctg_pos;
    

        
    // 去掉标题；
    getline(is, suftmp);
    string chrom_now;
    while(!is.eof()){
        is>>chrom;
        is>>recount;
        is>>length;
        {
            int pos = chrom.find('.', 0);
            chrom_now = chrom.substr(0, pos);
        }
        if (chrom_now.compare(chrom_group) == 0){
            LenMap.insert(pair<string ,int>(chrom, 1));
            s_set.push_back(chrom);
        }
        
    }
    
    is.clear();
    is.close();
    return 0;
}


int read_new_score(string input_file, vector<vector<float>> &dist_maxtrix){
    ifstream is(input_file, ios::in);
    if (!is.good()){
        cerr<<"the file is error!"<<endl;
        exit(-1);
    }
    string suftmp;
    string chrom_1;
    string chrom_2;
    int ctg_pos;
    int ctg_pos_invers;
    float score;

    while(!is.eof()){
        is>>chrom_1;
        is>>chrom_2;
        is>>score;
        {
            int pos = chrom_1.find("ctg", 0);
            ctg_pos=atoi((chrom_1.substr(pos+3, chrom_1.size()-pos-4).c_str()));
            pos = chrom_2.find("ctg", pos);
            ctg_pos_invers = atoi((chrom_2.substr(pos+3, chrom_2.size()-pos-3).c_str()));
            
        }
        // cout<<ctg_pos<<'\t'<<ctg_pos_invers<<endl;
        dist_maxtrix[ctg_pos][ctg_pos_invers] = score;
        dist_maxtrix[ctg_pos_invers][ctg_pos] = score;
    }

    is.close();
    return 0;
}


// 重载
int read_new_score(string input_file, vector<vector<float>> &dist_maxtrix, map<pair<int, int>, pair<char, char>> &OrientMap){
    ifstream is(input_file, ios::in);
    if (!is.good()){
        cerr<<"the "<<input_file<<" is error!"<<endl;
        exit(-1);
    }
    string suftmp;
    string chrom_1;
    string chrom_2;
    int ctg_pos;
    int ctg_pos_invers;
    float score;

    while(!is.eof()){
        is>>chrom_1;
        is>>chrom_2;
        is>>score;
        {
            int pos = chrom_1.find("ctg", 0);
            ctg_pos=atoi((chrom_1.substr(pos+3, chrom_1.size()-pos-4).c_str()));
            pos = chrom_2.find("ctg", pos);
            ctg_pos_invers = atoi((chrom_2.substr(pos+3, chrom_2.size()-pos-3).c_str()));
            
        }
        //cout<<ctg_pos<<'\t'<<ctg_pos_invers<<endl;
        dist_maxtrix[ctg_pos][ctg_pos_invers] = score;
        dist_maxtrix[ctg_pos_invers][ctg_pos] = score;
        OrientMap.insert(pair<pair<int, int>, pair<char, char> >({{ctg_pos, ctg_pos_invers}, {chrom_1[chrom_1.size()-1], chrom_2[chrom_2.size()-1]}}));
        OrientMap.insert(pair<pair<int, int>, pair<char, char> >({{ctg_pos_invers, ctg_pos}, {chrom_2[chrom_2.size()-1], chrom_1[chrom_1.size()-1]}}));

    }

    is.close();
    return 0;
}





int get_compare_iter(string start_chrom, map<float, vector<string> > &ResMap, vector<vector<float>> &dist_maxtrix, string chrom_group){

    int len = 1;
    int ctg_len = TMP_LEN;


    vector<map<string, float> > all_sort_arr(ctg_len);

    int ctg_pos = 0;
    int ctg_pos_invers = 0;


    map<string, int> s_v_set=LenMap;
    
    int ctg_first;
    {
        int pos = start_chrom.find("ctg", 0);
        ctg_first=atoi((start_chrom.substr(pos+3, start_chrom.size()-pos-3).c_str()));
    }
    vector<float> dist = dist_maxtrix[ctg_first];


    vector<string> v_set;
    string ctg_name;
    v_set.push_back(start_chrom);
    s_v_set[start_chrom] = 0;
    //float max = *max_element(dist.begin(), dist.end());
    //int max_index = max_element(dist.begin(), dist.end()) - dist.begin();
    //ctg_name = "Chr1.ctg"+to_string(max_index);

    

    float max;
    int max_index;
    ctg_name = start_chrom;

    while(v_set.size() != LenMap.size()){
        if (s_v_set[ctg_name] == 0){
            float max_s_v = 0;
            int max_s_v_index = 0;
            string pos_name;
            string pos_index;
            for(int i=0; i < dist.size(); i++){
                if (dist[i] == 0) continue;
                pos_index=chrom_group+".ctg"+to_string(i);
                if (max_s_v == 0){
                    if ((s_v_set[pos_index] == 1)){
                        max_s_v = dist[i];
                        max_s_v_index = i;
                        pos_name=chrom_group+".ctg"+to_string(i);
                    }
                }else{
                    if ((max_s_v < dist[i]) && (s_v_set[pos_index] == 1)){
                        max_s_v = dist[i];
                        max_s_v_index = i;
                        pos_name=chrom_group+".ctg"+to_string(i);
                    }
                }
                
            }
            max = max_s_v;
            max_index = max_s_v_index;
            ctg_name = pos_name;
        }
        v_set.push_back(ctg_name);

        
        s_v_set[ctg_name] = 0;
        string suftmp;
        for(int i=0; i < dist_maxtrix[max_index].size(); i++){
            //cout<<it->first<<'\t'<<it->second<<endl;
            if (dist_maxtrix[max_index][i] == 0) continue;

            suftmp = chrom_group+".ctg"+to_string(i);

            if (dist[i] == 0){
                dist[i] = max+dist_maxtrix[max_index][i];
            }else{
                if (dist[i] < (max + dist_maxtrix[max_index][i])){
                    dist[i] = max + dist_maxtrix[max_index][i];
                }
            }

        }

    }


    float S_total = 0;
    for(int i=0; i < v_set.size()-1; i++){
        int pair_i = 0;
        int pair_j = 0;

        {
            int pos = v_set[i].find("ctg", 0);
            int nex_pos = v_set[i+1].find("ctg", 0);
            if ((v_set[i] == "") || (v_set[i+1] == "")) continue;
            pair_i=atoi((v_set[i].substr(pos+3, v_set[i].size()-pos-3).c_str()));
            pair_j = atoi((v_set[i+1].substr(pos+3, v_set[i+1].size()-pos-3).c_str()));
            S_total += dist_maxtrix[pair_i][pair_j];
        }
    }
    ResMap.insert(pair<float, vector<string> >({S_total, v_set}));
    


    return 0;
}

int get_compare_iter(string start_chrom, map<float, vector<string> > &ResMap, vector<vector<float>> &dist_maxtrix, vector<vector<string>> &Orient_maxtrix, string chrom_group){

    int len = 1;
    int ctg_len = TMP_LEN;


    vector<map<string, float> > all_sort_arr(ctg_len);



    int ctg_pos = 0;
    int ctg_pos_invers = 0;

    map<string, int> s_v_set=LenMap;
    
    int ctg_first;
    {
        int pos = start_chrom.find("ctg", 0);
        ctg_first=atoi((start_chrom.substr(pos+3, start_chrom.size()-pos-3).c_str()));
    }
    vector<float> dist = dist_maxtrix[ctg_first];

    vector<string> orient_set;
    vector<string> v_set;
    string ctg_name;
    v_set.push_back(start_chrom);
    s_v_set[start_chrom] = 0;
    //float max = *max_element(dist.begin(), dist.end());
    //int max_index = max_element(dist.begin(), dist.end()) - dist.begin();
    //ctg_name = "Chr1.ctg"+to_string(max_index);

    

    float max;
    int max_index;
    ctg_name = start_chrom;

    while(v_set.size() != LenMap.size()){
        if (s_v_set[ctg_name] == 0){
            float max_s_v = 0;
            int max_s_v_index = 0;
            string pos_name;
            string pos_index;
            for(int i=0; i < dist.size(); i++){
                if (dist[i] == 0) continue;
                pos_index=chrom_group+".ctg"+to_string(i);
                if (max_s_v == 0){
                    if ((s_v_set[pos_index] == 1)){
                        max_s_v = dist[i];
                        max_s_v_index = i;
                        pos_name=chrom_group+".ctg"+to_string(i);
                    }
                }else{
                    if ((max_s_v < dist[i]) && (s_v_set[pos_index] == 1)){
                        max_s_v = dist[i];
                        max_s_v_index = i;
                        pos_name=chrom_group+".ctg"+to_string(i);
                    }
                }
                
            }
            max = max_s_v;
            max_index = max_s_v_index;
            ctg_name = pos_name;
        }
        v_set.push_back(ctg_name);

        
        s_v_set[ctg_name] = 0;
        string suftmp;
        for(int i=0; i < dist_maxtrix[max_index].size(); i++){
            //cout<<it->first<<'\t'<<it->second<<endl;
            if (dist_maxtrix[max_index][i] == 0) continue;

            suftmp = chrom_group+".ctg"+to_string(i);

            if (dist[i] == 0){
                dist[i] = max+dist_maxtrix[max_index][i];
            }else{
                if (dist[i] < (max + dist_maxtrix[max_index][i])){
                    dist[i] = max + dist_maxtrix[max_index][i];
                }
            }

        }

    }


    float S_total = 0;
    for(int i=0; i < v_set.size()-1; i++){
        int pair_i = 0;
        int pair_j = 0;

        {
            int pos = v_set[i].find("ctg", 0);
            int nex_pos = v_set[i+1].find("ctg", 0);
            if ((v_set[i] == "") || (v_set[i+1] == "")) continue;
            pair_i=atoi((v_set[i].substr(pos+3, v_set[i].size()-pos-3).c_str()));
            pair_j = atoi((v_set[i+1].substr(pos+3, v_set[i+1].size()-pos-3).c_str()));
            S_total += dist_maxtrix[pair_i][pair_j];
        }
    }
    ResMap.insert(pair<float, vector<string> >({S_total, v_set}));
    


    return 0;
}







int get_all_compare(vector<vector<float>> &dist_maxtrix, string chrom_group){

    int ctg_pos;

    string start_chrom;
    vector<string> Res;
    float S_res=0;
    map<float, vector<string>>ResMap;
    int pos = 0;

    string tmp_file = "/data0/stu_wangfang/tmp2/tmp_guancha.file";
    ofstream ot(tmp_file, ios::out);


    for(auto it=LenMap.rbegin(); it != LenMap.rend(); it++){
        start_chrom = it->first;
        get_compare_iter(start_chrom, ResMap, dist_maxtrix, chrom_group);
        for(auto itt=ResMap.begin(); itt != ResMap.end(); itt++){
            if (S_res == 0){
                S_res = itt->first;
                Res = itt->second;
            }else{
                if (S_res < itt->first){
                    S_res = itt->first;
                    Res = itt->second;
                }
            }
        }
          
        ResMap.clear();
    }
    
    
    // 优化步骤，可删除不影响
    // 效果不理想：因为插入的时候破坏了最大的顺序
    // vector<string> repeat_circle;
    // map<string, pair<string, string>> pairMap;
    // for(int i=0; i < Res.size()-1; i++){
    //     int max_index;
    //     int pair_i = 0;
    //     int pair_j = 0;
    //     {
    //          if ((Res[i] == "") || (Res[i+1] == "")) continue;
    //         int pos = Res[i].find("ctg", 0);
    //         pair_i=atoi((Res[i].substr(pos+3, Res[i].size()-pos-3).c_str()));
    //         pos = Res[i+1].find("ctg", 0);
    //         pair_j=atoi((Res[i+1].substr(pos+3, Res[i+1].size()-pos-3).c_str()));
    //     }
        
    //     if (dist_maxtrix[pair_i][pair_j] != 0){
    //         repeat_circle.push_back(Res[i]);
    //     }else{
    //         //cout<<Res[i+1]<<endl;
    //         max_index = max_element(dist_maxtrix[pair_j].begin(), dist_maxtrix[pair_j].end()) - dist_maxtrix[pair_j].begin();
    //         //cout<<max_index<<endl;
    //         {
    //             string ctg_max = "Chr1.ctg"+to_string(max_index);
    //             auto ctg_max_index = find(Res.begin(), Res.end(), ctg_max) - Res.begin();
    //             {
    //                 int pos = Res[ctg_max_index-1].find("ctg", 0);
    //                 int pos_pre=atoi((Res[ctg_max_index-1].substr(pos+3, Res[ctg_max_index-1].size()-pos-3).c_str()));
    //                 pos = Res[ctg_max_index+1].find("ctg", 0);
    //                 int pos_next=atoi((Res[ctg_max_index+1].substr(pos+3, Res[ctg_max_index+1].size()-pos-3).c_str()));
    //                 if (dist_maxtrix[max_index][pos_pre] < dist_maxtrix[max_index][pos_next]){
    //                     pairMap.insert(pair<string, pair<string, string>>({Res[ctg_max_index], {Res[i+1], Res[ctg_max_index]}}));
    //                 }else{
    //                     pairMap.insert(pair<string, pair<string, string>>({Res[ctg_max_index], {Res[ctg_max_index], Res[i+1]}}));
    //                 }
    //             }
    //         }
    //     }

    // }

    // vector<string>result;
    // map<string, int> s_v_set;
    // for(int i=0; i < s_set.size(); i++){
    //     s_v_set.insert(pair<string, int>(s_set[i], 1));
    // }
    // for(int i=0; i < Res.size(); i++){
    //     if ((pairMap.find(Res[i]) == pairMap.end()) && (s_v_set[Res[i]] == 1)){
    //         s_v_set[Res[i]] = 0;
    //         result.push_back(Res[i]);
    //     }else{
    //         if (s_v_set[Res[i]] == 1){
    //             result.push_back(pairMap[Res[i]].first);
    //             s_v_set[pairMap[Res[i]].first] = 0;
    //             result.push_back(pairMap[Res[i]].second);
    //             s_v_set[pairMap[Res[i]].second] = 0;
    //         }
    //     }
    // }


    // // 优化步骤到这里结束


    for(int i=0; i < Res.size(); i++){
        cout<<Res[i]<<endl;
    }
    
    cout<<S_res<<endl;

    ot.close();
    return 0;

}

// 重载
int get_all_compare(vector<vector<float>> &dist_maxtrix, map<pair<int, int>, pair<char, char> > &OrientMap, string chrom_group){

    int ctg_pos;

    string start_chrom;
    vector<string> Res;
    float S_res=0;
    map<float, vector<string>>ResMap;
    int pos = 0;

    string tmp_file = "/data0/stu_wangfang/tmp2/tmp_guancha.file";
    ofstream ot(tmp_file, ios::out);


    for(auto it=LenMap.rbegin(); it != LenMap.rend(); it++){
        start_chrom = it->first;
        get_compare_iter(start_chrom, ResMap, dist_maxtrix, chrom_group);
        for(auto itt=ResMap.begin(); itt != ResMap.end(); itt++){
            if (S_res == 0){
                S_res = itt->first;
                Res = itt->second;
            }else{
                if (S_res < itt->first){
                    S_res = itt->first;
                    Res = itt->second;
                }
            }
        }
          
        ResMap.clear();
    }

    // 纠正方向
    //vector<string> Result;
    vector<string> Cor_Result;
    //vector<int> Plus;
    vector<char> Plus_flag; 
    string chrom_1;
    string chrom_2;
    string chrom_1_new;
    string chrom_2_new;

    for(int _pos = 0; _pos < Res.size()-1; _pos++){
        int site = Res[_pos].find("ctg", 0);
        
        int ctg_first=atoi((Res[_pos].substr(site+3, Res[_pos].size()-site-3).c_str()));
        int rsite = Res[pos+1].find("ctg", 0);
        int rctg_first=atoi((Res[_pos+1].substr(rsite+3, Res[_pos+1].size()-rsite-3).c_str()));

        if (OrientMap.find({ctg_first, rctg_first}) == OrientMap.end()){
            chrom_1 = Res[_pos] + '+';
            chrom_2 = Res[_pos+1] + '-';
        }else{
            chrom_1 = Res[_pos]+OrientMap[{ctg_first, rctg_first}].first;
            chrom_2 = Res[_pos+1]+OrientMap[{ctg_first, rctg_first}].second;
        }

        Cor_Result.push_back(chrom_1);

        if (_pos == (Res.size()+2)){
            Cor_Result.push_back(chrom_2);
        }

    }

    string input_file = "/data0/stu_wangfang/tmp2/so_score.tour";
    ofstream ott(input_file, ios::out);
    for(int i=0; i < Cor_Result.size(); i++){
        ott<<Cor_Result[i]<<" ";
    }
    ott<<endl;
    ott.close();
    
    // //vector<string> Result;
    // vector<string> Cor_Result;
    // //vector<int> Plus;
    // vector<char> Plus_flag;
    // string chrom_1;
    // string chrom_2;
    // for(int i=0; i < Res.size()-1; i=i+2){
    //     {
    //         int pos = Res[i].find("ctg", 0);
    //         int ctg_first=atoi((Res[i].substr(pos+3, Res[i].size()-pos-3).c_str()));
    //         int rpos = Res[i+1].find("ctg", 0);
    //         int rctg_first=atoi((Res[i+1].substr(pos+3, Res[i+1].size()-pos-3).c_str()));
    //         chrom_1 = Res[i]+OrientMap[{ctg_first, rctg_first}].first;
    //         chrom_2 = Res[i+1]+OrientMap[{ctg_first, rctg_first}].second;
    //         cout<<OrientMap[{ctg_first, rctg_first}].first<<'\t'<<OrientMap[{ctg_first, rctg_first}].second<<endl;

    //         if (i==0){
    //             switch (chrom_1[chrom_1.size()-1])
    //             {
    //                 case '+':
    //                     Plus_flag.push_back('+');
    //                     break;
    //                 case '-':
    //                     Plus_flag.push_back('-');
    //                     break;
    //                 default:
    //                     break;
    //             }

    //             switch (chrom_2[chrom_2.size()-1])
    //             {
    //                 case '+':
    //                     Plus_flag.push_back('+');
    //                     break;
    //                 case '-':
    //                     Plus_flag.push_back('-');
    //                     break;
    //                 default:
    //                     break;
    //             }
    //             Cor_Result.push_back(chrom_1);
    //             Cor_Result.push_back(chrom_2);

    //         }else{
    //             if (chrom_1[chrom_1.size()-1] != Plus_flag[0]){
    //                 chrom_1[chrom_1.size()-1] = Plus_flag[0];
    //             }
    //             if (chrom_2[chrom_2.size()-1] != Plus_flag[1]){
    //                 chrom_2[chrom_2.size()-1] = Plus_flag[1];
    //             }
    //             Cor_Result.push_back(chrom_1);
    //             Cor_Result.push_back(chrom_2);
    //         }
    //     }
    // }

    /*
    // 纠正方向
    //vector<string> Result;
    vector<string> Cor_Result;
    //vector<int> Plus;
    vector<char> Plus_flag; 
    string chrom_1;
    string chrom_2;
    string chrom_1_new;
    string chrom_2_new;

    // 忘记考虑两个点直接不相连的情况了， sos
    int site = Res[0].find("ctg", 0);
    int ctg_first=atoi((Res[0].substr(site+3, Res[0].size()-site-3).c_str()));
    int rsite = Res[1].find("ctg", 0);
    int rctg_first=atoi((Res[1].substr(rsite+3, Res[1].size()-rsite-3).c_str()));
    chrom_1 = Res[0]+OrientMap[{ctg_first, rctg_first}].first;
    chrom_2 = Res[1]+OrientMap[{ctg_first, rctg_first}].second;

    Cor_Result.push_back(chrom_1);
    Cor_Result.push_back(chrom_2);
    chrom_1_new = chrom_1;
    chrom_2_new = chrom_2;
    
    int _pos = 1;
    while (_pos < Res.size()-1){
        
        site = Res[_pos].find("ctg", 0);
        int ctg_first=atoi((Res[_pos].substr(site+3, Res[_pos].size()-site-3).c_str()));
        int rsite = Res[pos+1].find("ctg", 0);
        int rctg_first=atoi((Res[_pos+1].substr(rsite+3, Res[_pos+1].size()-rsite-3).c_str()));

        chrom_1 = Res[_pos]+OrientMap[{ctg_first, rctg_first}].first;
        chrom_2 = Res[_pos+1]+OrientMap[{ctg_first, rctg_first}].second;
        chrom_2_new = chrom_2;

        if (chrom_1[chrom_1.size()-1] != chrom_1_new[chrom_1_new.size()-1]){
            if (chrom_2[chrom_2.size()-1] == '+'){
                chrom_2_new[chrom_2_new.size()-1] = '-';
            }else{
                chrom_2_new[chrom_2_new.size()-1] = '+';
            }
        }
        Cor_Result.push_back(chrom_2_new);
        _pos++;
        chrom_1_new=chrom_2_new;
    }

    */
            
        

    for(int i=0; i < Res.size(); i++){
        cout<<Res[i]<<endl;
    }
    
    cout<<S_res<<endl;


    ot.close();
    return 0;

}





int main(){
    //string input_clm = "/data2/wangyb/0.ALLPhase/2.pore-c/CEN/2.contig1/minimap2/test_optimize/AT_pore_c_1.clm";
    string input_file = "/data2/wangyb/0.ALLPhase/2.pore-c/CEN/2.contig1/minimap2/test_optimize/AT_pore_c_1.counts_GATC.group1.txt";

    string input_new_score = "/data0/stu_wangfang/tmp2/so_score.txt";

    vector<vector<float>> dist_maxtrix(TMP_LEN);
    for(int i=0; i < TMP_LEN; i++){
        dist_maxtrix[i].resize(TMP_LEN);
    }
    map<pair<int, int>, pair<char, char>> OrientMap;    

    string chrom_group = "Chr1";
    read_new_score(input_new_score, dist_maxtrix, OrientMap);
    get_distance(input_file, chrom_group);
    get_all_compare(dist_maxtrix, OrientMap, chrom_group);

    // string pair_score_file = "/data2/wangyb/0.ALLPhase/2.pore-c/CEN/2.contig1/minimap2/test_optimize/new.score2.txt";
    // string new_max_file = "/data0/stu_wangfang/tmp2/new_score2.txt";
    // get_max_pair(pair_score_file, new_max_file);

    return 0;
}