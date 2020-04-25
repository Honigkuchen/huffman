#ifndef HUFFMAN_TABLE_H
#define HUFFMAN_TABLE_H

// STL includes
#include <vector>

// Project includes
#include "types.h"

namespace huffman
{
/*!
 * \brief The Table class holds the Huffman codes for each individual symbol.
 *
 * The template parameter S is the symbol type and must be comparable to be inserted into a map.
 */
template<typename S, typename F = long long, typename = std::enable_if_t<std::is_signed_v<F>>>
class Table
{
public:
	using SymbolType = S;
	using FrequencyType = F;
	struct TableEntry
	{
		TableEntry(const SymbolType& s, const FrequencyType& f, const BinaryNumber& n) :
			symbol(s),
			frequency(f),
			binary_number(n) {}

		SymbolType symbol;
		FrequencyType frequency;
		BinaryNumber binary_number;

		void PrintOn(std::ostream& o) const
		{
			o << "'" << symbol << "': " << frequency << " -> ";
			binary_number.PrintOn(o);
			o << std::endl;
		}
	};
	using EntryType = TableEntry;
public:
    /*!
     * \brief Table Creates an empty table.
     */
	constexpr Table() noexcept = default;
	/*!
	 *
	 */
	constexpr Table(std::initializer_list<EntryType>&& v) noexcept :
		symbols_(std::move(v)) {}
    /*!
     * \brief Table Creates a table with the passed data.
     * \param data The new data of the table
     */
    template<typename T>
    explicit constexpr Table(const T& data) :
        symbols_(data) {}
    /*!
     * \brief Table Creates a table with the passed data.
	 *
	 * The data is moved from the argument into the internal data structure thus being destroyed.
	 *
     * \param data The new data of the table
     */
    template<typename T>
    explicit constexpr Table(T&& data) noexcept :
        symbols_(std::move(data)) {}
    /*!
     * \brief operator [] Lookup operator for the table.
	 *
     * \param symbol The symbol to get the Huffman code from

     * \return The Huffman binary number of the symbol
     */
    [[nodiscard]] constexpr BinaryNumber operator[](const SymbolType& symbol) const
    {
        if(auto iter = std::find_if(symbols_.begin(), symbols_.end(), [&symbol](const EntryType& entry)
		{ return entry.symbol == symbol; }); iter != symbols_.end())
            return iter->binary_number;
        else return BinaryNumber();
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
	/*!
	 *
	 */
	constexpr void Append(const Table<S, F>& other)
	{
		symbols_.insert(symbols_.end(), other.symbols_.begin(), other.symbols_.end());
	}
private:
    std::vector<EntryType> symbols_;
};
}

#endif // HUFFMAN_TABLE_H
