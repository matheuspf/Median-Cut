#include <bits/stdc++.h>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#define USING_OPENCV

#include "../Wrapper/Operations.h"
#include "MedianCut.h"
#include "../Benchmark.h"



/* g++ -I"C:\opencv\build\include" -I"C:\MinGW\include" -L"C:\opencv\mingw64\lib" main.cpp C:\MinGW\lib\libboost_thread.a C:\MinGW\lib\libboost_system.a  -lopencv_core249 -lopencv_highgui249 -o main.exe -std=c++14 -O3 */



using namespace std;
using namespace cv;

struct Foo
{
    template <typename T, typename U>
    inline auto operator () (const wrp::Operations<cv::Vec<T, 3>>& v, const wrp::Operations<cv::Vec<U, 3>>& u)
    //inline auto operator () (const cv::Vec<T, 3>& v, const cv::Vec<U, 3>& u)
    {
        return std::abs(v[0] - u[0]) + std::abs(v[1] - u[1]) + std::abs(v[2] - u[2]);
    }
};




int main()
{
    Mat_<Vec3b> img(imread("..//img1.jpg"));

    vector<wrp::Operations<Vec3b>> v;
    //vector<Vec3b> v;

    for(int i = 0; i < img.rows; ++i)
        for(int j = 0; j < img.cols; ++j)
            v.push_back(img(i, j));


    //cout << Benchmark<milli>([&]{ mc::MedianCutVariance<wrp::Operations<Vec3b>, 3, 32, Foo>()(v); }) << endl;
    //cout << Benchmark<milli>([&]{ mc::MedianCutVariance<Vec3b, 3, 32, Foo>()(v); }) << endl;


    auto res = mc::MedianCutDefault<wrp::Operations<Vec3b>, 3, 32, dist::Euclidean>()(v);

    for(int i = 0; i < get<0>(res).size(); ++i)
    {
        for(auto x : get<1>(res)[i])
        {
            for(int j = 0; j < 3; ++j)
            {
                //cout << int(get<0>(res)[i][j]) << " ";
                img(x / img.cols, x % img.cols)[j] = get<0>(res)[i][j];
            }
            //cout << endl;
        }
    }

    namedWindow("w");
    imshow("w", img);
    waitKey(0);


    return 0;
}


