#ifndef HUFFMAN_H
#define HUFFMAN_H

#ifndef __cplusplus
#error This library can only be used in C++.
#else
#if !(__cplusplus >= 201703L || _MSVC_LANG >= 201703L)
#error C++17 is needed for this library.
#endif
#endif

// STL includes
#include <vector>
#include <deque>
#include <map>
#include <memory>
#include <algorithm>

// Project includes
#include "binary_number.h"
#include "table.h"
#include "graph.h"

namespace huffman
{
/*!
 * \brief TraverseTree This function creates the binary codes by traversing the Huffman graph.
 *
 * The traversal is done recursively. If the graph is too large (like really really large...) this may lead to a stack overflow, pay attention!
 *
 * \param node The current node to traverse from
 * \param prefix The current prefix to append the current and all child nodes to
 *
 * \return A map of symbols and their corresponding Huffman binary representation
 */
template<typename S, typename F = long long, typename = std::enable_if_t<std::is_signed_v<F>>>
[[nodiscard]] Table<S, F> TraverseTree(const std::unique_ptr<Node_<F>>& node, const BinaryNumber& prefix)
{
	using SymbolType = S;
	using FrequencyType = F;

	using NodeType = Node_<FrequencyType>;
	using InternalNodeType = InternalNode_<FrequencyType>;
	using LeafNodeType = LeafNode_<SymbolType, FrequencyType>;

    if(InternalNodeType* in = dynamic_cast<InternalNodeType*>(node.get()); in != nullptr)
    {
		Table<S, F> codes;
        BinaryNumber left_prefix = prefix;
        left_prefix.AppendBack(BinaryDigit::ZERO);
		codes.Append(TraverseTree<SymbolType>(in->left, left_prefix));

        BinaryNumber right_prefix = prefix;
        right_prefix.AppendBack(BinaryDigit::ONE);
		codes.Append(TraverseTree<SymbolType>(in->right, right_prefix));
        return std::move(codes);
    }
	else if (LeafNodeType* ln = dynamic_cast<LeafNodeType*>(node.get()); ln != nullptr)
	{
		BinaryNumber internal_prefix;
		if (prefix.empty()) internal_prefix.AppendBack(BinaryDigit::ZERO);
		else internal_prefix = prefix;
		const Table<S, F> codes = { {ln->symbol, ln->frequency, internal_prefix} };
		return std::move(codes);
	}
    else
    {
        // Should never happen!
        throw std::runtime_error("ERROR: This should never happen, unknown Node_ Pointer detected. Aborting.");
    }
}

/*!
 * \brief Encode This function creates Huffman codes for the passed symbols.

 * \param symbols The symbols to get the Huffman codes for

 * \return A table of codes for each individual symbol
 */
template<typename S, typename F = long long, typename = std::enable_if_t<std::is_signed_v<F>>>
[[nodiscard]] Table<S> Encode(const std::vector<S>& symbols)
{
	using SymbolType = S;
	using FrequencyType = F;

	using NodeType = Node_<FrequencyType>;
	using NodeTypePtr = std::unique_ptr<NodeType>;
	using InternalNodeType = InternalNode_<FrequencyType>;
	using LeafNodeType = LeafNode_<SymbolType, FrequencyType>;

    std::map<SymbolType, FrequencyType> frequencies;

    for(const SymbolType& symbol : symbols)
    {
        if(auto iter = frequencies.find(symbol); iter != frequencies.end())
            ++iter->second;
        else
            frequencies.insert(std::make_pair(symbol, 1));
    }

    std::deque<NodeTypePtr> nodes;

	for (const auto&[symbol, frequency] : frequencies)
		nodes.push_back(std::make_unique<LeafNodeType>(symbol, frequency));
	
    constexpr auto NodeComparator = [](const NodeTypePtr& lhs, const NodeTypePtr& rhs)
    {
        return lhs->frequency < rhs->frequency;
    };

    std::sort(nodes.begin(), nodes.end(), NodeComparator);

    while(nodes.size() > 1)
    {
		NodeTypePtr right_child = std::move(nodes.front());
        nodes.pop_front();
		NodeTypePtr left_child = std::move(nodes.front());
        nodes.pop_front();
        nodes.push_front(std::make_unique<InternalNodeType>(right_child, left_child));
        std::sort(nodes.begin(), nodes.end(), NodeComparator);
    }

	NodeTypePtr root = std::move(nodes.front());

    return TraverseTree<SymbolType>(root, {});
}

}

#endif // HUFFMAN_H
