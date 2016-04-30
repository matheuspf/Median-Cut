#ifndef MEDIAN_CUT_H
#define MEDIAN_CUT_H

#include <type_traits>
#include <tuple>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

#include "../Distances/Distances.h"



namespace mc
{

namespace impl
{

template <typename, std::size_t, std::size_t, class>
struct DefaultCut;

template <typename, std::size_t, std::size_t, class>
struct VarianceCut;


//--------------------------------------------------------------------------------------


// Returns at compile time if the type T can be multiplied by type U

template <typename T, typename U>
constexpr bool isMultipliable (void*, std::decay_t<decltype(std::declval<T>() * std::declval<U>())>* = nullptr) { return true; }

constexpr bool isMultipliable (...) { return false; }



// With the helpers bellow, one can easily do 'TypePrec_t<T, U>' without worrying if
// type T can be multiplied with a type U
//
// If you have, for example, a std::array<char, N>, which can not be multiplied with a double,
// the traits returns the std::array<char, N>
//
// But if you create a wrapper, you can define the result of a multiplication with a double
// to be std::array<double, N>, for example

template <typename T, typename U, bool>
struct TypePrecImpl
{
    using Type = decltype(std::declval<T>() * std::declval<U>());
};

template <typename T, typename U>
struct TypePrecImpl<T, U, false>
{
    using Type = T;
};


template <typename T, typename U>
struct TypePrec : TypePrecImpl<std::decay_t<T>, std::decay_t<U>, true> {};


template <typename T, typename U>
using TypePrec_t = typename TypePrec<T, U>::Type;


}




//========================================================================================================



template < typename T, std::size_t Dim, std::size_t K, class Distance = dist::Euclidean,
           template <typename, std::size_t, std::size_t, class> class Cutter = impl::DefaultCut
         >
class MedianCut : public Cutter<T, Dim, K, Distance>
{
public:

    using Base = Cutter<T, Dim, K, Distance>;

    using Base::Value;
    using Base::Cut;


    MedianCut (const Distance& dist = Distance()) : Base(dist) {}


    inline auto operator () (const std::vector<T>& data)
    {
        return this->operator()(data, std::vector<double>(data.size(), 1));
    }

    auto operator () (const std::vector<T>& data, const std::vector<double>& weights)
    {
        std::array<std::vector<int>, K> points = { std::vector<int>(data.size()) };

        std::iota(points[0].begin(), points[0].end(), 0);



        for(int k = 1; k < K; ++k)
        {
            double val = std::numeric_limits<double>::min();
            int pos = 0;

            for(int i = 0; i < k; ++i)
            {
                double temp = Value(data, points[i], weights);

                if(val < temp)
                {
                    val = temp;
                    pos = i;
                }
            }

            points[k] = Cut(data, points[pos], weights);
        }


        std::array<impl::TypePrec_t<T, double>, K> centroids;


        for(int i = 0; i < K; ++i)
        {
            impl::TypePrec_t<T, double> aux{};

            for(int x : points[i])
            {
                for(int j = 0; j < Dim; ++j)
                    aux[j] += data[x][j];
            }

            for(int j = 0; j < Dim; ++j)
                centroids[i][j] = aux[j] / points[i].size();
        }


        return std::make_tuple(centroids, points);
    }
};


template <typename T, std::size_t Dim, std::size_t K, class Distances>
using MedianCutDefault = MedianCut<T, Dim, K, Distances, impl::DefaultCut>;

template <typename T, std::size_t Dim, std::size_t K, class Distances>
using MedianCutVariance = MedianCut<T, Dim, K, Distances, impl::VarianceCut>;



namespace impl
{

template <typename T, std::size_t Dim, std::size_t K, class Distance>
class DefaultCut
{
public:

    DefaultCut (const Distance& dist = Distance()) : dist(dist) {}



    auto Value (const std::vector<T>& data, const std::vector<int>& points, const std::vector<double>& weights)
    {
        return std::get<0>(AuxFunc(data, points, weights));
    }


    auto Cut (const std::vector<T>& data, std::vector<int>& points, const std::vector<double>& weights)
    {
        int dim = std::get<1>(AuxFunc(data, points, weights));

        std::nth_element(points.begin(), points.begin() + points.size() / 2, points.end(),
                         [&data, dim](int x, int y){ return data[x][dim] < data[y][dim]; });


        std::vector<int> r;

        std::move(points.begin() + points.size() / 2, points.end(), std::inserter(r, r.end()));

        points.resize(points.size() / 2);


        return r;
    }

private:

    Distance dist;



protected:

    auto AuxFunc (const std::vector<T>& data, const std::vector<int>& points, const std::vector<double>&)
    {
        T lower, upper;

        for(int i = 0; i < Dim; ++i)
        {
            lower[i] = std::numeric_limits<std::decay_t<decltype(lower[i])>>::max();
            upper[i] = std::numeric_limits<std::decay_t<decltype(upper[i])>>::min();
        }

        for(auto x : points)
        {
            for(int i = 0; i < Dim; ++i)
            {
                lower[i] = std::min(data[x][i], lower[i]);
                upper[i] = std::max(data[x][i], upper[i]);
            }
        }


        auto distance = upper[0] - lower[0];
        int pos = 0;

        for(int i = 1; i < Dim; ++i)
        {
            auto aux = upper[i] - lower[i];

            if(distance < aux)
            {
                distance = aux;
                pos = i;
            }
        }


        return std::make_tuple(std::move(distance), pos);
    }
};


template <typename T, std::size_t Dim, std::size_t K, class Distance>
class VarianceCut : public DefaultCut<T, Dim, K, Distance>
{
public:

    VarianceCut (const Distance& dist = Distance()) : dist(dist) {}


    auto Value (const std::vector<T>& data, const std::vector<int>& points, const std::vector<double>& weights)
    {
        impl::TypePrec_t<T, double> mean{};
        double variance = 0.0, div = 0.0;

        for(auto p : points)
        {
            for(int i = 0; i < Dim; ++i)
               mean[i] += data[p][i];  // * weights[p];

            //div += weights[p];
        }

        for(int i = 0; i < Dim; ++i)
            mean[i] /= points.size(); //div;


        for(int p : points)
            variance += dist(data[p], mean); // * weights[p];


        return variance;
    }

private:

    Distance dist;
};


}   // namespace impl


}  // namespace mc


////=============================================================================================



#endif // MEDIAN_CUT_H
