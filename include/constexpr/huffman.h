#pragma once

// STL includes
#include <iostream>
#include <limits>

// Project includes
#include "table.h"
#include "binary_number.h"

// Project includes
#include "types.h"

namespace huffman
{
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