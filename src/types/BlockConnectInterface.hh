#ifndef __BLOCKCONNECTINTERFACE_HH
#define __BLOCKCONNECTINTERFACE_HH

#include "types/Boundary.hh"

#include <upcxx/upcxx.hpp>

template <typename T>
struct BlockConnectInterface {
	Boundary boundary;
	int size;
	int stride;
	int startIndex;
	T pointerH;
	T pointerB;
	T pointerHu;
	T pointerHv;
};

#endif // __BLOCKCONNECTINTERFACE_HH