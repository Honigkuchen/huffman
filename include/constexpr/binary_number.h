#pragma once

#define ENABLE_IOSTREAM 1

// C includes
#include <cmath>

// STL includes
#include <algorithm>
#include <memory>

// Project includes
#include "types.h"

namespace huffman
{
/*!
 * \brief enum class BinaryDigit This enum represents the two binary digits zero and one.
 */
enum class BinaryDigit
{
    ZERO = 0,
    ONE = 1,
    NONE = 2
};

[[nodiscard]] constexpr auto to_underlying(BinaryDigit d) noexcept
{
    return static_cast<std::underlying_type_t<BinaryDigit>>(d);
}
/*!
 * \brief PrintBinaryDigit This function prints a binary digit on the specified output stream.
 *
 * \param d The digit to print
 * \param o The output stream to write to
 */
constexpr void PrintBinaryDigit(const BinaryDigit& d, std::ostream& o)
{
#if ENABLE_IOSTREAM
    switch (d)
    {
    case BinaryDigit::ZERO:
        o << 0;
        break;
    case BinaryDigit::ONE:
        o << 1;
        break;
    case BinaryDigit::NONE:
        o << '-';
        break;
    }
#endif
}
/*!
 * \brief PrintByte This function prints a byte in binary representation on the specified output stream.
 *
 * \param b The byte to print
 * \param o The output stream to write to
 */
constexpr void PrintByte(const Byte& b, std::ostream& o)
{
#if ENABLE_IOSTREAM
    for (auto i = 7; i >= 0; --i)
    {
        auto mask = 1 << i;
        if ((b & mask) != 0)
            PrintBinaryDigit(BinaryDigit::ONE, o);
        else
            PrintBinaryDigit(BinaryDigit::ZERO, o);
    }
#endif
}
/*!
 * \brief class BinaryNumber 
 */
 template<std::size_t container_size>
class BinaryNumber
{
public:
    void PrintOn(std::ostream& o) const
    {
#if ENABLE_IOSTREAM
        for (const auto& digit : digits)
            PrintBinaryDigit(digit, o);
#endif
    }
    constexpr void AppendBack(const BinaryDigit& d)
    {
        //static_assert(index < container_size - 1, "Can not append binary digit to binary value!");
        if(index < container_size - 1)
        {
            digits[index] = d;
            ++index;
        }
    }
    constexpr void AppendBack(const unsigned short& d)
    {
        if(d == 0 || d == 1)
            AppendBack(BinaryDigit(d));
    }
    template<typename T>
    constexpr void AppendBack(std::initializer_list<T>&& new_digits)
    {
        for(const auto& d : new_digits)
            AppendBack(d);
    }
    [[nodiscard]] constexpr auto begin() noexcept
    {
        return digits.begin();
    }
    [[nodiscard]] constexpr auto begin() const noexcept
    {
        return digits.begin();
    }
    [[nodiscard]] constexpr auto end() noexcept
    {
        return digits.end();
    }
    [[nodiscard]] constexpr auto end() const noexcept
    {
        return digits.end();
    }
    [[nodiscard]] constexpr auto empty() const noexcept
    {
        return digits.empty();
    }
    constexpr void fill(const BinaryDigit& v)
    {
        digits.fill(v);
    }
    [[nodiscard]] constexpr auto size() const noexcept
    {
        return digits.size();
    }
    [[nodiscard]] const auto& operator[](const std::size_t& pos) const
    {
        return digits[pos];
    }
    [[nodiscard]] const auto& at(const std::size_t& pos) const
    {
        return digits.at(pos);
    }
    bool operator==(const BinaryNumber& other) const noexcept
    {
        if (digits.size() != other.digits.size()) return false;
        for (auto this_iter = digits.rbegin(), other_iter = other.digits.rbegin();
            this_iter != digits.rend() && other_iter != other.digits.rend();
            ++this_iter, ++other_iter)
        {
            if (*this_iter != *other_iter) return false;
        }
        return true;
    }
    void swap(BinaryNumber& rhs) noexcept
    {
        digits.swap(rhs.digits);
    }
    friend void swap(BinaryNumber& first, BinaryNumber& second) noexcept
    {
        first.swap(second);
    }
private:
    std::size_t index = 0;
    std::array<BinaryDigit, container_size> digits;
};


}
