#define ENABLE_IOSTREAM 1

// C includes
#include <cmath>

// STL includes
#include <algorithm>
#include <variant>
#include <memory>
#include <vector>
#if ENABLE_IOSTREAM
#include <iostream>
#endif
#include <deque>

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

using Byte = char;

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
    std::vector<Byte> ToByteRepresentation() const noexcept
    {
        constexpr auto BitsPerByte = 8;
        std::deque<BinaryDigit> local_copy = digits;

        auto remainder = BitsPerByte - (local_copy.size() % BitsPerByte);
        std::deque<BinaryDigit> expanded_binary_number(remainder);
        for (auto i = 0; i < remainder; ++i)
            expanded_binary_number[i] = BinaryDigit::NONE;

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
    std::size_t index = 0;
    std::array<BinaryDigit, container_size> digits;
};

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
     * \brief operator [] Lookup operator for the table.
	 *
     * \param symbol The symbol to get the Huffman code from
     * \return The Huffman binary number of the symbol
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

template<typename S = char, typename F = long long>
class Node
{
public:
    void PrintOn(std::ostream& o) const
    {
#if ENABLE_IOSTREAM
        o << symbol << ": " << frequency << " -> " << (isValid ? "isValid" : "isNotValid") << "\n";
#endif
    }
    void PrintRecursiveOn(std::ostream& o) const
    {
#if ENABLE_IOSTREAM
        if(isValid)
        {
            PrintOn(o);
            if(left != nullptr)
                left->PrintRecursiveOn(o);
            if(right != nullptr)
                right->PrintRecursiveOn(o);
        }
#endif
    }
    F frequency = std::numeric_limits<F>::min();
    S symbol = '\0';
    Node* left = nullptr;
    Node* right = nullptr;
    bool isValid = false;
    bool isLeafNode = true;
};

using NodeType = Node<>;

constexpr auto count_valid_nodes(auto begin, auto end)
{
    auto result = 0;
    for(auto iter = begin; iter != end; ++iter)
    {
        if(iter->isValid)
            ++result;
    }
    return result;
}

constexpr auto get_first_two_valid_nodes(auto begin, auto end)
{
    bool first_picked = false;
    std::pair<NodeType*, NodeType*> result;
    bool both_found = false;
    for(auto iter = begin; iter != end; ++iter)
    {
        if(iter->isValid)
            if(!first_picked)
            {
                result.first = &(*iter);
                first_picked = true;
            }
            else
            {
                result.second = &(*iter);
                both_found = true;
            }
        if(both_found)
            break;
    }
    return result;
}

template<std::size_t count, std::size_t return_count = 2 * count, std::size_t double_return_count = 2 * return_count>
[[nodiscard]] constexpr const Table<count> TraverseTree(const NodeType& node, BinaryNumber<count> prefix)
{
    using SymbolType = char;
	using FrequencyType = long long;
    using BinaryNumberType = decltype(prefix);

    if(!node.isLeafNode)
    {
        BinaryNumberType left_prefix = prefix;
        left_prefix.AppendBack(BinaryDigit::ZERO);
        const Table<count> left_traversed_tree = TraverseTree(*node.left, left_prefix);

        BinaryNumberType right_prefix = prefix;
        right_prefix.AppendBack(BinaryDigit::ONE);
        const Table<count> right_traversed_tree = TraverseTree(*node.right, right_prefix);

        auto codes = Merge(left_traversed_tree, right_traversed_tree);

        return codes;
    }
	else
	{
		BinaryNumberType internal_prefix;
		if (prefix.empty())
            internal_prefix.AppendBack(BinaryDigit::ZERO);
		else
            internal_prefix = prefix;
		Table<count> codes;
        codes.add_entry(node.symbol, node.frequency, internal_prefix);
		return codes;
	}
}

//template<std::size_t count, typename S = char, typename F = long long>
[[nodiscard]] constexpr const auto Encode(const auto symbols)
{
    static_assert(std::is_literal_type<decltype(symbols)>::value, "symbols is not a literal!");

    using SymbolType = char;
	using FrequencyType = long long;
    const std::size_t count = symbols.size();
    std::array<std::pair<SymbolType, FrequencyType>, count> frequencies;

    frequencies.fill({'\0', std::numeric_limits<FrequencyType>::min()});
    
    for(auto counter = 0; const SymbolType& symbol : symbols)
    {
        if(auto iter = std::find_if(frequencies.begin(), frequencies.end(), [&symbol](const auto& f)
        {
            return f.first == symbol;
        }); iter != frequencies.end())
            ++iter->second;
        else
        {
            frequencies[counter].first = symbol;
            frequencies[counter].second = 1;
            ++counter;
        }
    }

    struct NodeTree
    {
        std::array<NodeType, 2*count> nodes;
        std::array<NodeType, 2*count> swap_nodes;
        constexpr auto& root() const
        {
            return nodes[0];
        }
    };

    NodeTree tree;

	for (auto counter = 0; const auto&[symbol, frequency] : frequencies)
    {
        if(frequency != std::numeric_limits<FrequencyType>::min())
        {
            auto& node = tree.nodes[counter];
            node.frequency = frequency;
            node.symbol = symbol;
            node.isValid = true;
            ++counter;
        }
    }

    constexpr auto NodeComparator = [](const NodeType& lhs, const NodeType& rhs)
    {
        return lhs.frequency > rhs.frequency && lhs.isValid;
    };

    std::sort(tree.nodes.begin(), tree.nodes.end(), NodeComparator); 

    auto counter = 0;
    auto insert_counter = count_valid_nodes(tree.nodes.begin(), tree.nodes.end());
    const auto amount_of_different_symbols = insert_counter;

    auto print_arrays = [&tree]()
    {
#if ENABLE_IOSTREAM
        constexpr bool print = true;
        if constexpr(print)
        {
            std::cout << "Nodes:\n";
            for(const auto& node : tree.nodes)
                node.PrintOn(std::cout);
            
            std::cout << "Swap Nodes:\n";
            for(const auto& node : tree.swap_nodes)
                node.PrintOn(std::cout);
            std::cout << "________________" << std::endl;
        }
#endif
    }; 

    print_arrays();

    while(count_valid_nodes(tree.nodes.begin(), tree.nodes.end()) > 1)
    {
        auto valid_nodes = get_first_two_valid_nodes(tree.nodes.begin(), tree.nodes.end());
        std::swap(tree.swap_nodes[counter], *valid_nodes.first);
        std::swap(tree.swap_nodes[counter + 1], *valid_nodes.second);
        valid_nodes.first->isValid = false;
        valid_nodes.second->isValid = false;

        auto& inserted_node = tree.nodes[insert_counter];
        inserted_node.isValid = true;
        inserted_node.isLeafNode = false;
        inserted_node.left = &tree.swap_nodes[counter];
        inserted_node.right = &tree.swap_nodes[counter + 1];
        inserted_node.frequency = inserted_node.left->frequency + inserted_node.right->frequency;

        std::sort(tree.nodes.begin(), tree.nodes.end(), NodeComparator);
        counter += 2;
        ++insert_counter;

        print_arrays();
    }
    std::sort(tree.nodes.begin(), tree.nodes.end(), NodeComparator);
    print_arrays();
    return TraverseTree(tree.root(), BinaryNumber<count>());
}

}

int main()
{
    constexpr std::array data {'a', 'b', 'a', 'c'};

    //static constinit auto table = huffman::Encode(data);
    auto table = huffman::Encode(data);
#if ENABLE_IOSTREAM
    table.PrintOn(std::cout);
#endif

    return int(table[1][1]);
}