#pragma once

#// Project includes
#include "binary_number.h"
#include "types.h"


namespace huffman 
{
template<typename T, std::size_t N1> 
[[nodiscard]] constexpr auto join_arrays(const std::array<T, N1> first, const std::array<T, N1> second)
{
    std::array<T, N1> c;

    auto counter = 0;

    for(const auto& element : first)
    {
        if(element.frequency != std::numeric_limits<long long>::max())
        {
            c[counter] = element;
            ++counter;
        }
    }

    for(const auto& element : second)
    {
        if(element.frequency != std::numeric_limits<long long>::max())
        {
            c[counter] = element;
            ++counter;
        }
    }

    return c;
}

/*!
 * \brief The Table class holds the Huffman codes for each individual symbol.
 *
 * The template parameter S is the symbol type and must be comparable to be inserted into a map.
 */
template<std::size_t container_size, typename S = char, typename F = long long>
class Table
{
public:
	using SymbolType = S;
	using FrequencyType = F;

    template<std::size_t binary_number_digits>
	struct TableEntry
	{
        constexpr TableEntry() noexcept :
            symbol('\0'),
            frequency(std::numeric_limits<F>::min())
        {
            binary_number.fill(BinaryDigit::NONE);
        }
		constexpr TableEntry(const SymbolType& s, const FrequencyType& f, const BinaryNumber<binary_number_digits>& n) :
			symbol(s),
			frequency(f),
			binary_number(n) {}

		SymbolType symbol = '\0';
		FrequencyType frequency = std::numeric_limits<F>::min();
		BinaryNumber<binary_number_digits> binary_number;

		void PrintOn(std::ostream& o) const
		{
#if ENABLE_IOSTREAM
			o << "'" << symbol << "': " << frequency << " -> ";
			binary_number.PrintOn(o);
			o << std::endl;
#endif
		}
	};
	using EntryType = TableEntry<container_size>;
    using BinaryNumberType = BinaryNumber<container_size>;
public:
    /*!
     * \brief operator [] Lookup operator for the
using Byte = char;
ymbol
     */
    [[nodiscard]] constexpr BinaryNumberType operator[](const SymbolType& symbol) const
    {
        if(auto iter = std::find_if(symbols_.begin(), symbols_.end(), [&symbol](const EntryType& entry)
		{ return entry.symbol == symbol; }); iter != symbols_.end())
            return iter->binary_number;
        else return BinaryNumberType();
    }
	/*!
	 * \brief PrintOn Prints the entire table on a specified output stream.
	 *
	 * \param o The output stream to write to
	 */
	constexpr void PrintOn(std::ostream& o) const noexcept
	{
		for (const auto& s : symbols_)
			s.PrintOn(o);
	}
    
    constexpr void add_entry(const auto& s, const auto& f, const auto& n)
    {
        symbols_[elements_index_].symbol = s;
        symbols_[elements_index_].frequency = f;
        symbols_[elements_index_].binary_number = n;
        ++elements_index_;
    }

    [[nodiscard]] constexpr std::size_t size() const noexcept
    {
        return container_size;
    }

    std::size_t elements_index_ = 0;
    std::array<EntryType, container_size> symbols_;
};

template<std::size_t N1>
[[nodiscard]] constexpr const Table<N1> Merge(const Table<N1> first, const Table<N1> second)
{
    Table<N1> result;
    result.symbols_ = join_arrays(first.symbols_, second.symbols_);
    return result;
}

}