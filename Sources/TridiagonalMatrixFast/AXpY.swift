//
//  AXpY.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

/// Computes `y ← a·x + y` in-place.
///
/// - Parameters:
///   - a: Scalar multiplier.
///   - x: Input vector.
///   - y: Vector to be updated.
@usableFromInline @inline(__always)
func aXpY_Inplace<T: ScalarField>(a: T, x: [T], y: inout [T]) {
	precondition(x.count == y.count, "Vector size mismatch")
	let n = Int32(x.count)
	if n == 0 { return }
	T.axpy(n, a, x, 1, &y, 1)
}

/// Computes `a·x + y` and returns the result.
///
/// - Parameters:
///   - a: Scalar multiplier.
///   - x: Input vector.
///   - y: Input vector.
/// - Returns: The sum `a·x + y`.
@usableFromInline @inline(__always)
func aXpY<T: ScalarField>(a: T, x: [T], y: [T]) -> [T] {
	var out = y
	aXpY_Inplace(a: a, x: x, y: &out)
	return out
}

// MARK: - Internal AXpY Implementations

/// Real scalar optimized implementation of `y ← A·x + y` for tridiagonal `A`.
@usableFromInline @inline(__always)
func AXpY_<T: RealScalar>(
	_ A: TridiagonalMatrix<T>, _ x: [T], _ y: inout [T]
) -> [T] {
	let n = x.count
	precondition(n == A.size)
	precondition(y.count == n)
	
	if n == 0 { return [] }
	
	if n == 1 {
		y[0] += A.diagonal[0] * x[0]
		return y
	}
	
	x.withUnsafeBufferPointer { xPtr in
		y.withUnsafeMutableBufferPointer { yPtr in
			// y += d * x
			A.diagonal.withUnsafeBufferPointer { dPtr in
				T.vma(dPtr.baseAddress!, 1,
					  xPtr.baseAddress!, 1,
					  yPtr.baseAddress!, 1,
					  n)
			}
			// y[0...n-2] += upper * x[1...n-1]
			A.upper.withUnsafeBufferPointer { uPtr in
				T.vma(uPtr.baseAddress!, 1,
					  xPtr.baseAddress! + 1, 1,
					  yPtr.baseAddress!, 1,
					  n - 1)
			}
			// y[1...n-1] += lower * x[0...n-2]
			A.lower.withUnsafeBufferPointer { lPtr in
				T.vma(lPtr.baseAddress!, 1,
					  xPtr.baseAddress!, 1,
					  yPtr.baseAddress! + 1, 1,
					  n - 1)
			}
		}
	}
	return y
}

/// Helper: complex multiply-add for a band: `y += band * x`.
@usableFromInline @inline(__always)
func complexVMA<T: RealScalar>(
	_ band: [Complex<T>],
	_ x: UnsafePointer<T>,
	_ y: CMutablePtr<T>,
	_ temp: CMutablePtr<T>,
	_ count: Int
) {
	band.withUnsafeBufferPointer { bandPtr in
		bandPtr.baseAddress!.withMemoryRebound(
			to: T.self,
			capacity: 2 * count
		) { b in
			T.vma(b, 2, x, 2, y, 2, count)          // y.real += b.real * x.real
			T.vmul(b + 1, 2, x + 1, 2, temp, 1, count)  // temp = b.imag * x.imag
			T.vsub(temp, 1, y, 2, y, 2, count)      // y.real = y.real - temp
			T.vma(b, 2, x + 1, 2, y + 1, 2, count)  // y.imag += b.real * x.imag
			T.vma(b + 1, 2, x, 2, y + 1, 2, count)  // y.imag += b.imag * x.real
		}
	}
}

/// Complex scalar implementation of `y ← A·x + y` for tridiagonal `A`.
@usableFromInline @inline(__always)
func AXpY_<T: RealScalar>(
	_ A: TridiagonalMatrix<Complex<T>>,
	_ x: [Complex<T>],
	_ y: inout [Complex<T>]
) -> [Complex<T>] {
	let n = x.count
	precondition(n == A.size)
	precondition(y.count == n)
	
	if n == 0 { return y }
	
	if n == 1 {
		y[0] = y[0] + A.diagonal[0] * x[0]
		return y
	}
	
	let tempSize = n
	
	if tempSize <= 1024 {
		withUnsafeTemporaryAllocation(of: T.self, capacity: tempSize) { temp in
			computeAXpYWithTemps(A, x: x, y: &y, temp: temp.baseAddress!)
		}
	} else {
		let temp = UnsafeMutableBufferPointer<T>.allocate(capacity: tempSize)
		defer { temp.deallocate() }
		computeAXpYWithTemps(A, x: x, y: &y, temp: temp.baseAddress!)
	}
	return y
}

/// Complex AXpY using a shared temporary buffer for diagonal and off-diagonals.
@inlinable func computeAXpYWithTemps<T: RealScalar>(
	_ A: TridiagonalMatrix<Complex<T>>,
	x: [Complex<T>],
	y: inout [Complex<T>],
	temp: UnsafeMutablePointer<T>
) {
	let n = x.count
	
	x.withUnsafeBufferPointer { xPtr in
		y.withUnsafeMutableBufferPointer { yPtr in
			xPtr.baseAddress!.withMemoryRebound(to: T.self, capacity: 2 * n) { x in
				yPtr.baseAddress!.withMemoryRebound(to: T.self, capacity: 2 * n) { y in
					complexVMA(A.diagonal, x, y, temp, n)
					complexVMA(A.upper, x + 2, y, temp, n - 1)
					complexVMA(A.lower, x, y + 2, temp, n - 1)
				}
			}
		}
	}
}
