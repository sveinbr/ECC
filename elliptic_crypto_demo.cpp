struct integer // Define as a separate type to spot automatic type conversions.
{
    typedef int value_type;
    value_type i;
    integer operator%(value_type p) const noexcept { return integer{ (i%p) + (i < 0 ? p : 0) }; } // Always 0 <= result < p
    integer operator*(integer other) const noexcept { return integer{ i*other.i }; }
    //integer operator+(integer other) const noexcept { return integer{ i + other.i }; }
    //integer operator-() const noexcept { return integer{ -i }; }
};
constexpr integer operator"" _i(unsigned long long i) { return integer{ static_cast<integer::value_type>(i) }; }

template <integer::value_type p>
integer InvertModuloP(integer x) noexcept
{
    auto ExtendedEuclidianAlgorithm = [](integer::value_type a, integer::value_type b)
    {
        struct { integer::value_type r, s, t; } new_{ b,0,1 }, old_{ a,1,0 };
        while (new_.r != 0)
        {
            integer::value_type quotient = old_.r / new_.r;
            decltype(old_) temp_ = old_;
            old_ = new_;
            new_ = { temp_.r - quotient*old_.r, temp_.s - quotient*old_.s, temp_.t - quotient*old_.t }; // new = temp - quotient*old
        }
        return old_;
    };

    auto divisor = ExtendedEuclidianAlgorithm(x.i % p, p);
    return integer{ divisor.r == 1 ? divisor.s : 0 } % p;
}



template <integer::value_type p>
struct FP_element
{
    integer value;
    bool operator==(FP_element other) const noexcept { return value.i == other.value.i; }
    FP_element operator-() const noexcept { return{ integer{ value.i==0 ? 0 : p - value.i } }; } // Optimized "(-value)%p"
    FP_element operator+(FP_element other) const noexcept  { return FP_element{ value.i + other.value.i >= p ? value.i + other.value.i - p : value.i + other.value.i }; }  // Optimized "(value+other.value)%p"
    FP_element operator-(FP_element other) const noexcept { return FP_element{ value.i >= other.value.i ? value.i - other.value.i : p - (other.value.i - value.i)}; } // Optimized "(value-other.value)%p"
    FP_element operator*(FP_element other) const noexcept { return FP_element{ (other.value*value) % p }; }
    FP_element operator/(FP_element divisor) const noexcept { return FP_element{ (value * InvertModuloP<p>(divisor.value)) % p }; }
    friend FP_element<p> operator*(integer k, FP_element<p> x) noexcept { return FP_element<p>{ (k*x.value) % p }; }
};



template <integer::value_type p, integer::value_type a, integer::value_type b>
struct Point
{
    FP_element<p> x, y;
    Point operator-() const noexcept { return Point{ x, -y }; }
    Point operator-(Point other) const noexcept { return *this + (-other); }
    bool operator==(Point other) const noexcept { return x == other.x && y == other.y; }
    friend Point<p, a, b> operator+(Point<p, a, b> P, Point<p, a, b> Q) noexcept
    {
        const Point<p, a, b> ZERO{ 0,0 };
        if (P == ZERO) return Q;
        if (Q == ZERO) return P;
        if (P.y == -Q.y) return{ 0,0 };
        FP_element<p> m = P == Q ?
            (integer{ 3 } *P.x * P.x + FP_element<p>{a}) / (integer{ 2 } *P.y)
            : (P.y - Q.y) / (P.x - Q.x);
        FP_element<p> x = m * m - P.x - Q.x;
        FP_element<p> y = P.y + m*(x - P.x);
        return{ x, -y };
    }
    friend Point<p, a, b> operator*(integer k, Point<p, a, b> P) noexcept
    {
        Point Q = P;
        while (--k.value.i > 0) // O(k) for simplicity, can be improved to O(log(k))
            Q = Q + P;
        return Q;
    }
};

// Defining a specific curve of the class y^2 = x^3 + ax + b in the space F_p
template <integer::value_type FieldOrder, integer::value_type A, integer::value_type B>
struct EllipticCurve
{
    static const integer::value_type a = A, b = B, p = FieldOrder;
    typedef ::Point<p, a, b> Point;

    bool contains(Point P) const noexcept
    {
		auto y_2 = P.y*P.y;
		auto x_3 = P.x*P.x*P.x;
		auto ax = FP_element<p>{ a }*P.x;
		auto b_ = FP_element<p>{ b };
		auto rhs = x_3 + ax + b_;
        return P.y*P.y == P.x*P.x*P.x + FP_element<p>{ a }*P.x + FP_element<p>{ b };
    }

    integer::value_type order(Point P) const noexcept
    {
        if (!contains(P))
            return 0;
        const Point ZERO{ 0,0 };
        Point nP = P;
        integer::value_type n{ 1 };
        while (!(nP == ZERO))
        {
            ++n;
            nP = nP + P;
        }
        return n;
    }

    integer::value_type order() const noexcept // brute force. Consider Schoof's algorithm.
    {
        integer::value_type out{ 1 }; // Include infinity
        for (integer::value_type x = 0; x < p; ++x)
            for (integer::value_type y = 0; y < p; ++y)
                if (contains(Point{ x,y }))
                    out++;
        return out;
    }
};


#include <iostream>
template <integer::value_type p> void Require_that_mutiplication_and_division_modulo_p_are_consistent()
{
    int errors = 0;
    for (integer::value_type x = 1; x < p; x++)
    {
        for (integer::value_type y = 1; y < p; y++)
        {
            FP_element<p> X{ x }, Y{ y };
            auto X2 = (X*Y) / Y, X3 = (X / Y) * Y;
            if (X.value.i != X2.value.i || X.value.i != X3.value.i)
            {
                std::cout << "FAIL: " << __FUNCTION__ << ", X=" << x << ", Y=" << y << std::endl;
                errors++;
            }
        }
    }
    if (errors == 0)
        std::cout << "OK:   " << __FUNCTION__ << std::endl;
}


void Require_that_EC_7_6_3_has_order_7_and_base_point_3_0_with_order_7()
{
    typedef EllipticCurve<7, 6, 3> EC;
    EC C;
	if (C.order() == 6)
		std::cout << "OK:   " << __FUNCTION__ << std::endl;
	else
		std::cout << "FAIL: " << __FUNCTION__ << "EC [p=" << C.p << ", a=" << C.a << ", b=" << C.b << "] has order " << C.order() << " (expected 6)" << std::endl;

	for (integer::value_type y = EC::p - 1; y >= 0; --y)
    {
		for (integer::value_type x = 0; x < EC::p; ++x)
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
    Require_that_mutiplication_and_division_modulo_p_are_consistent<7>();
    Require_that_EC_7_6_3_has_order_7_and_base_point_3_0_with_order_7();
    return 0;
}
