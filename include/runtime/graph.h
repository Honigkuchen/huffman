#ifndef HUFFMAN_GRAPH_H
#define HUFFMAN_GRAPH_H

// STL includes
#include <memory>

namespace huffman
{
	/*!
	 * \brief The Node_ class represents a node in the Huffman graph.
	 *
	 * Each node has a frequency for the represented symbol or for the intermediate node.
	 * The frequency is represented by the template parameter T. The type must
	 * meet the std::is_signed requirement.
	 *
	 * @see LeafNode_
	 * @see InternalNode_
	 */
	template<typename F = long long, typename = std::enable_if_t<std::is_signed_v<F>>>
	class Node_
	{
	public:
		/*! \brief The type of the frequency. */
		using FrequencyType = F;
	public:
		/*!
		 * \brief Node_ Creates a Node_ object with a specified frequency for the represented symbol.
		 *
		 * If the symbol is lower than zero a runtime error is thrown.
		 *
		 * \param frequency The frequency of the symbol
		 */
		explicit constexpr Node_(const FrequencyType& frequency) :
			frequency(frequency)
		{
			//if (frequency < 0)
			//	throw std::runtime_error("The frequency of a node must not be negative");
		}
		// Prevent slicing
		Node_(const Node_&) = delete;
		Node_& operator=(const Node_&) = delete;
		/*!
		 * \brief ~Node_ Default Node_ destructor.
		 */
		virtual ~Node_() = default;
	public:
		const FrequencyType frequency; //! frequency This is the frequency of the symbol or of the intermediate node in thje graph.
	};
	/*!
	 * \brief Default node type, default frequency type is 'long long'.
	 */
	using Node = Node_<>;
	/*!
	 * \brief The InternalNode_ class represents an internal node in the Huffman graph.
	 *
	 * An InternalNode_ does not represent a single symbol. The assigned frequency is the sum of the child nodes.
	 * The template parameter FrequencyType is 'long long' by default and must meet the std::is_signed requirement.
	 *
	 * @see Node_
	 * @see LeafNode_
	 */
	template<typename F = long long>
	class InternalNode_ : public Node_<F>
	{
	public:
		/*!
		 * \brief InternalNode_ Creates an internal node in the Huffman graph by setting the child nodes.
		 * \param left The left child node
		 * \param right The right child node
		 */
		constexpr InternalNode_(std::unique_ptr<Node_>& left, std::unique_ptr<Node_>& right) :
			Node(left->frequency + right->frequency),
			left(std::move(left)),
			right(std::move(right)) {}
	public:
		const std::unique_ptr<Node> left; //! left The left child node
		const std::unique_ptr<Node> right; //! right The right child node
	};
	/*!
	 * \brief Default internal node type, default frequency type is 'long long'.
	 */
	using InternalNode = InternalNode_<>;
	/*!
	 * \brief The LeafNode_ class represents the leafes of the Huffman graph.
	 *
	 * Each leaf node represents a single symbol. Each symbol occurs with a spcified frequency.
	 * The frequency type must meet the std::is_signed requirement.
	 *
	 * @see Node_
	 * @see InternalNode_
	 */
	template<typename S, typename F = long long>
	class LeafNode_ : public Node_<F>
	{
	public:
		/*! \brief The type of the symbol. */
		using SymbolType = S;
		/*! \brief The type of the frequency. */
		using FrequencyType = F;
	public:
		/*!
		 * \brief LeafNode_ Creates a leaf node for a given symbol and frequency.
		 * \param symbol The symbol represented by this node
		 * \param frequency The frequency of the symbol
		 */
		constexpr LeafNode_(const SymbolType& symbol, const FrequencyType& frequency) :
			Node_(frequency),
			symbol(symbol) {}
	public:
		const SymbolType symbol; //! symbol The symbol represented by this node
	};
	template<typename SymbolType>
	using LeafNode = LeafNode_<SymbolType>;
}

#endif // HUFFMAN_GRAPH_H
