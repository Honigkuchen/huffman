#include "huffman.h"

int main()
{
    constexpr std::array data = {'a', 'a', 'b', 'd'};
    //static constinit auto table = huffman::Encode(data);
    auto table = huffman::Encode(data);
    return 0;
}