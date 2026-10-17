#ifndef MONTGOMERY_HPP
#define MONTGOMERY_HPP

#include <cassert>
#include <cstdint>

struct Montgomery64 {
	using U = std::uint64_t;
	using W = __uint128_t;
	U mod, inverse = 1, r2, one;
	explicit Montgomery64(U n) : mod(n), r2(-W(n) % n) {
		for (int repeat_index = 0; repeat_index < 6; ++repeat_index) inverse *= 2 - n * inverse;
		assert(n > 1 && (n & 1));
		one = reduce(r2);
	}
	U reduce(W x) const {
		U a = x >> 64, b = (W(U(x) * inverse) * mod) >> 64;
		return a < b ? a - b + mod : a - b;
	}
	U in(U x) const { return reduce(W(x) * r2); }
	U mul(U a, U b) const { return reduce(W(a) * b); }
	U add(U a, U b) const {
		U c = a + b;
		return c < a || c >= mod ? c - mod : c;
	}
	U pow(U a, U e) const {
		U r = one;
		for (; e; e >>= 1, a = mul(a, a))
			if (e & 1) r = mul(r, a);
		return r;
	}
};

#endif
