#ifndef HUFFMAN_BINARY_NUMBER_H
#define HUFFMAN_BINARY_NUMBER_H

// STL includes
#include <queue>
#include <iostream>

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
		ONE = 1
	};
	/*!
	 * \brief PrintBinaryDigit This function prints a binary digit on the specified output stream.
	 *
	 * \param d The digit to print
	 * \param o The output stream to write to
	 */
	constexpr void PrintBinaryDigit(const BinaryDigit& d, std::ostream& o)
	{
		switch (d)
		{
		case BinaryDigit::ZERO:
			o << 0;
			break;
		case BinaryDigit::ONE:
			o << 1;
			break;
		}
	}
	/*!
	 * \brief PrintByte This function prints a byte in binary representation on the specified output stream.
	 *
	 * \param b The byte to print
	 * \param o The output stream to write to
	 */
	constexpr void PrintByte(const Byte& b, std::ostream& o)
	{
		for (auto i = 7; i >= 0; --i)
		{
			auto mask = 1 << i;
			if ((b & mask) != 0)
				PrintBinaryDigit(BinaryDigit::ONE, o);
			else
				PrintBinaryDigit(BinaryDigit::ZERO, o);
		}
	}
	/*!
	 * \brief class BinaryNumber 
	 */
	class BinaryNumber
	{
	public:
		BinaryNumber() = default;
		BinaryNumber(const BinaryNumber& other) = default;
		BinaryNumber& operator=(const BinaryNumber& other) = default;
		BinaryNumber(BinaryNumber&& other) noexcept = default;
		BinaryNumber& operator=(BinaryNumber&& other) noexcept = default;

		template<typename T>
		constexpr BinaryNumber(std::initializer_list<T>&& new_digits) noexcept
		{
			AppendBack(std::move(new_digits));
		}
		void PrintOn(std::ostream& o) const
		{
			for (const auto& digit : digits)
				PrintBinaryDigit(digit, o);
		}
		void AppendBack(const BinaryDigit& d)
		{
			digits.push_back(d);
		}
		void AppendBack(const unsigned short& d)
		{
			if(d == 0 || d == 1)
				digits.push_back(BinaryDigit(d));
		}
		template<typename T>
		constexpr void AppendBack(std::initializer_list<T>&& new_digits)
		{
			for(const auto& d : new_digits)
				AppendBack(d);
		}
		[[nodiscard]] auto begin() noexcept
		{
			return digits.begin();
		}
		[[nodiscard]] auto begin() const noexcept
		{
			return digits.begin();
		}
		[[nodiscard]] auto end() noexcept
		{
			return digits.end();
		}
		[[nodiscard]] auto end() const noexcept
		{
			return digits.end();
		}
		[[nodiscard]] auto empty() const noexcept
		{
			return digits.empty();
		}
		[[nodiscard]] auto size() const noexcept
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
		std::vector<Byte> ToByteRepresentation() const noexcept
		{
			constexpr auto BitsPerByte = 8;
			std::deque<BinaryDigit> local_copy = digits;

			auto remainder = BitsPerByte - (local_copy.size() % BitsPerByte);
			std::deque<BinaryDigit> expanded_binary_number(remainder);
			for (auto i = 0; i < remainder; ++i)
				expanded_binary_number[i] = BinaryDigit::ZERO;

			auto iterator = expanded_binary_number.begin() + remainder;

			expanded_binary_number.insert(iterator, local_copy.begin(), local_copy.end());

			std::vector<Byte> result;
			result.reserve(static_cast<std::size_t>(std::ceil(expanded_binary_number.size() / BitsPerByte)));

			for (auto i = 0; i < expanded_binary_number.size(); i += BitsPerByte)
			{
				Byte b = 0b00000000;
				for (auto j = BitsPerByte - 1; j >= 0; --j)
					if (expanded_binary_number[i + j] == BinaryDigit::ONE)
						b += static_cast<Byte>(std::pow(2, (BitsPerByte - 1 - j)));
					
				result.push_back(b);
			}
			return result;
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
		std::deque<BinaryDigit> digits;
	};
}

#endif // HUFFMAN_BINARY_NUMBER_H
