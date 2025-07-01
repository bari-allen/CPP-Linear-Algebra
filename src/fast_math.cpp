#include <cstdint>
#include <limits>
#include <bit>

namespace fast_math {

    double fast_sqrt(const double number) {
        static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 must be used!");

        if (number < 0) {return std::numeric_limits<double>::quiet_NaN();}

        constexpr double smallest_normal = std::numeric_limits<double>::min();
        if (number < smallest_normal) { return 0;}


        double d_number = number;
        double half = d_number * 0.5;
        //magic = 3/2 * 2^52 * (1023 - 0.043)
        constexpr uint64_t magic = 0x5FE339A0219FF02D;
        constexpr double three_halves = 1.5;

        //Evil bit hack
        uint64_t reinterpret = std::bit_cast<uint64_t>(d_number);
        reinterpret = magic - (reinterpret >> 1);
        d_number = std::bit_cast<double>(reinterpret);

        //Newton iterations to get a closer approximation
        d_number = d_number * (three_halves - (half * d_number * d_number));
        d_number = d_number * (three_halves - (half * d_number * d_number));
        d_number = d_number * (three_halves - (half * d_number * d_number));
        d_number = d_number * (three_halves - (half * d_number * d_number));

        //S * S ^ -1/2 = S ^ 1/2
        return d_number * number;
    }

    double fast_sqrt(const int number) {
        static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 must be used!");

        if (number < 0) {return std::numeric_limits<double>::quiet_NaN();}


        double d_number = static_cast<double>(number);
        double half = d_number * 0.5;
        //magic = 3/2 * 2^52 * (1023 - 0.043)
        constexpr uint64_t magic = 0x5FE339A0219FF02D;
        constexpr double three_halves = 1.5;

        //Evil bit hack
        uint64_t reinterpret = std::bit_cast<uint64_t>(d_number);
        reinterpret = magic - (reinterpret >> 1);
        d_number = std::bit_cast<double>(reinterpret);

        //Newton iterations to get a closer approximation
        d_number = d_number * (three_halves - (half * d_number * d_number));
        d_number = d_number * (three_halves - (half * d_number * d_number));
        d_number = d_number * (three_halves - (half * d_number * d_number));
        d_number = d_number * (three_halves - (half * d_number * d_number));

        //S * S ^ -1/2 = S ^ 1/2
        return d_number * static_cast<double>(number);
    }

    double abs(const double number) {
        uint64_t reinterpret = std::bit_cast<uint64_t>(number);
        constexpr uint64_t mask = 0x7FFFFFFFFFFFFFFF;
        reinterpret &= mask;

        double result = std::bit_cast<double>(reinterpret);
        return result;
    }

    uint32_t abs(const int32_t number) {
        int32_t mask = number >> 31;
        int32_t abs = (number ^ mask) - mask;

        return abs;
    }

}