

#include <iostream>
#include <fstream>
#include <iomanip>
#include "CycloidGear.h"

using namespace std;

int main()
{
    // 输出文件地址
    string file_name("E:\\md.txt");
    ofstream fs("E:\\md.txt");

    vector< vector<double>> md;

    double D0 = 177;
    double d0 = 7;
    double e = 1.7;
    int Z = 39;
    double dm = 8.7312;
    int n = 3;

    CycloidGear rv125N;

    rv125N.setZ(Z);
    rv125N.setrm(dm);
    rv125N.sete(e);

    double delta = 0.0001;
    double d125n = d0 -int(n / 2) * delta;     //   7.01 (-0.01 , +0.02)  :       7~7.03    :   0.00976 (161.995: 7 - 7.03
    double D125n = D0 -int(n / 2) * delta; // 162.00 (-0.005, +0.005) : 161.995~162.005 : + 0.00966 (7: 161.995-162.005)
                            //                                             = 0.01942 ((161.995, 7) - (162.005, 7.03))
    fs << 0 << " ";

    double di = d125n, d = d125n;
    double Di = D125n, D = D125n;

    for (int i = 0; i < n; i++, di += delta) {
        vector<double> mdR;
        mdR.push_back(di);
        for (int j = 0; j < 1; j++, Di += 0.001)
        {
            rv125N.setR(Di);
            rv125N.setr(di);
            auto bar = rv125N.find_bar();
            double md_R_r = rv125N.printMd(bar.first);
            mdR.push_back(md_R_r);

            printf("[%d, %d]bar: (%.5f, %f)\n", i, j, bar.first, bar.second);
        }
        md.push_back(mdR);
    }

    // 写入表头至文件
    for (int i = 0; i < md.size(); i++, D += delta)
    {
        fs << std::setprecision(10) << D << " ";
    }

    fs << "\n";

    
    // 写入数据至文件
    for (int i = 0; i < md.size(); ++i) {
        for (int j = 0; j < md[0].size(); ++j)
            fs << std::setprecision(10) << md[i][j] << " ";
        fs << "\n";
    }

    fs.close();

    return 0;
}
