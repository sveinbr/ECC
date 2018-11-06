template<typename Int>
Int modulus(Int i, Int p) noexcept { return (i%p) + (i < 0 ? p : 0); } // Always 0 <= modulus(i,p) < p, thus not the same as operator% for negative values

template<typename Int, Int p>
Int InvertModuloP(Int x) noexcept
{
    auto ExtendedEuclidianAlgorithm = [](Int a, Int b)
    {
        struct { Int r, s, t; } new_{ b,0,1 }, old_{ a,1,0 };
        while (new_.r != 0)
        {
            Int quotient = old_.r / new_.r;
            decltype(old_) temp_ = old_;
            old_ = new_;
            new_ = { temp_.r - quotient*old_.r, temp_.s - quotient*old_.s, temp_.t - quotient*old_.t }; // new = temp - quotient*old
        }
        return old_;
    };

    auto divisor = ExtendedEuclidianAlgorithm( modulus(x, p), p);
    return divisor.r == 1 ? modulus(divisor.s, p) : 0;
}


template<typename Int, Int p>
struct FP_element
{
    Int value;
    bool operator==(FP_element other) const noexcept { return value == other.value; }
	FP_element operator-() const noexcept { return FP_element{ value == 0 ? 0 : p - value }; } // Optimized "modulus(-value,p)"
    FP_element operator+(FP_element other) const noexcept  { return FP_element{ value + other.value >= p ? value + other.value - p : value + other.value }; }  // Optimized "modulus(value+other.value, p)"
    FP_element operator-(FP_element other) const noexcept { return FP_element{ value >= other.value ? value - other.value : p - (other.value - value)}; } // Optimized "modulus(value-other.value, p)"
    FP_element operator*(FP_element other) const noexcept { return FP_element{ modulus(other.value * value, p) }; }
    FP_element operator/(FP_element divisor) const noexcept {  return FP_element{ modulus(value * InvertModuloP<Int, p>(divisor.value), p) };}
    friend FP_element<Int,p> operator*(Int k, FP_element<Int,p> x) noexcept { return FP_element<Int,p>{ modulus(k*x.value, p) }; }
};



template <typename Int, Int p, Int a, Int b>
struct Point
{
    FP_element<Int, p> x, y;
    Point operator-() const noexcept { return Point{ x, -y }; }
    Point operator-(Point other) const noexcept { return *this + (-other); }
    bool operator==(Point other) const noexcept { return x == other.x && y == other.y; }
    friend Point<Int, p, a, b> operator+(Point<Int, p, a, b> P, Point<Int, p, a, b> Q) noexcept
    {
        const Point<Int, p, a, b> ZERO{ 0,0 };
        if (P == ZERO) return Q;
        if (Q == ZERO) return P;
        if (P.y == -Q.y) return{ 0,0 };
        FP_element<Int, p> m = P == Q ?
            (Int{ 3 } * P.x * P.x + FP_element<Int, p>{a}) / (Int{ 2 } * P.y)
            : (P.y - Q.y) / (P.x - Q.x);
        FP_element<Int, p> x = m * m - P.x - Q.x;
        FP_element<Int, p> y = P.y + m*(x - P.x);
        return{ x, -y };
    }
    friend Point<Int, p, a, b> operator*(Int k, Point<Int, p, a, b> P) noexcept
    {
        Point Q = P;
        while (--k.value.i > 0) // O(k) for simplicity, can be improved to O(log(k))
            Q = Q + P;
        return Q;
    }
};

// Defining a specific curve of the class y^2 = x^3 + ax + b in the space F_p
template <typename Int, Int FieldOrder, Int A, Int B>
struct EllipticCurve
{
    static const Int a = A, b = B, p = FieldOrder;
    typedef ::Point<Int, p, a, b> Point;
    typedef FP_element<Int, p> Element;

    bool contains(Point P) const noexcept
    {
        return P.y*P.y == P.x*P.x*P.x + a*P.x + Element{ b };
    }

    Int order(Point P) const noexcept
    {
        if (!contains(P))
            return 0;
        const Point ZERO{ 0,0 };
        Point nP = P;
        Int n{ 1 };
        while (!(nP == ZERO))
        {
            ++n;
            nP = nP + P;
        }
        return n;
    }

    Int order() const noexcept // brute force. Consider Schoof's algorithm.
    {
        Int out{ 1 }; // Include infinity
        for (Int x = 0; x < p; ++x)
            for (Int y = 0; y < p; ++y)
                if (contains(Point{ x,y }))
                    out++;
        return out;
    }
};


#include <iostream>
template <typename Int, Int p> void Require_that_mutiplication_and_division_modulo_p_are_consistent()
{
    int errors = 0;
    for (Int x = 1; x < p; x++)
    {
        for (Int y = 1; y < p; y++)
        {
            FP_element<Int,p> X{ x }, Y{ y };
            auto X2 = (X*Y) / Y, X3 = (X / Y) * Y;
            if (X.value != X2.value || X.value != X3.value)
            {
                std::cout << "FAIL: " << __FUNCTION__ << ", X=" << x << ", Y=" << y << std::endl;
                errors++;
            }
        }
    }
    if (errors == 0)
        std::cout << "OK:   " << __FUNCTION__ << std::endl;
}


template <typename Int> void Require_that_EC_7_6_3_has_order_7_and_base_point_3_0_with_order_7()
{
    typedef EllipticCurve<Int, 7, 6, 3> EC;
    EC C;
    if (C.order() == 6)
        std::cout << "OK:   " << __FUNCTION__ << std::endl;
    else
        std::cout << "FAIL: " << __FUNCTION__ << "EC [p=" << C.p << ", a=" << C.a << ", b=" << C.b << "] has order " << C.order() << " (expected 6)" << std::endl;

    for (Int y = EC::p - 1; y >= 0; --y)
    {
        for (Int x = 0; x < EC::p; ++x)
        {
            if (C.contains({ x,y }))
                std::cout << C.order({ x,y }) << "   ";
            else
                std::cout << ".   ";
        }
        std::cout << std::endl;
    }
}

int main()
{
    Require_that_mutiplication_and_division_modulo_p_are_consistent<int,7>();
    Require_that_EC_7_6_3_has_order_7_and_base_point_3_0_with_order_7<int>();
    Require_that_mutiplication_and_division_modulo_p_are_consistent<long, 7>();
    Require_that_EC_7_6_3_has_order_7_and_base_point_3_0_with_order_7<long>();
    return 0;
}
