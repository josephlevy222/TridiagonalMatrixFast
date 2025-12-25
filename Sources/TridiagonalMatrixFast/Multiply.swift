//
//  Multiply.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

/// Computes the matrix–vector product `A·x` for real scalars.
///
/// - Parameters:
///   - A: Tridiagonal matrix.
///   - x: Column vector.
/// - Returns: The product `A·x`.
@inlinable
public func multiply_<T: RealScalar>(
	_ A: TridiagonalMatrix<T>, _ x: [T]
) -> [T] {
	precondition(x.count == A.size, "Invalid column vector size")
	let n = x.count
	if n == 0 { return [] }
	if n == 1 { return [A.diagonal[0] * x[0]] }
	
	var result = [T](unsafeUninitializedCapacity: n) { buffer, initializedCount in
		initializedCount = n
	}
	
	x.withUnsafeBufferPointer { xPtr in
		result.withUnsafeMutableBufferPointer { resultPtr in
			let xBase = xPtr.baseAddress!
			let r = resultPtr.baseAddress!
			A.diagonal.withUnsafeBufferPointer { dPtr in
				T.vmul(dPtr.baseAddress!, 1,
					   xBase, 1,
					   r, 1, n)
			}
			A.upper.withUnsafeBufferPointer { uPtr in
				T.vma(uPtr.baseAddress!, 1,
					  xBase + 1, 1,
					  r, 1, n - 1)
			}
			A.lower.withUnsafeBufferPointer { lPtr in
				T.vma(lPtr.baseAddress!, 1,
					  xBase, 1,
					  r + 1, 1, n - 1)
			}
		}
	}
	
	return result
}

/// Computes the matrix–vector product `A·x` for complex scalars.
///
/// - Parameters:
///   - A: Complex tridiagonal matrix.
///   - x: Complex column vector.
/// - Returns: The product `A·x`.
@inlinable
public func multiply_<T: RealScalar>(
	_ A: TridiagonalMatrix<Complex<T>>,
	_ x: [Complex<T>]
) -> [Complex<T>] {
	precondition(x.count == A.size, "Invalid column vector size")
	let n = x.count
	if n == 0 { return [] }
	if n == 1 { return [A.diagonal[0] * x[0]] }
	
	var result = [Complex<T>](unsafeUninitializedCapacity: n) { buffer, initializedCount in
		initializedCount = n
	}
	
	var temp = [T](unsafeUninitializedCapacity: n) { buffer, initializedCount in
		initializedCount = n
	}
	
	x.withUnsafeBufferPointer { xPtr in
		result.withUnsafeMutableBufferPointer { resultPtr in
			xPtr.baseAddress!.withMemoryRebound(to: T.self, capacity: 2 * n) { x in
				resultPtr.baseAddress!.withMemoryRebound(to: T.self, capacity: 2 * n) { y in
					A.diagonal.withUnsafeBufferPointer { dPtr in
						dPtr.baseAddress!.withMemoryRebound(
							to: T.self,
							capacity: 2 * n
						) { d in
							T.vmul(d, 1, x, 1, y, 1, 2 * n)  // y.real = d.real*x.real ; y.imag = d.imag*x.imag
							T.vsub(y + 1, 2, y, 2, y, 2, n)  // y.real = y.real - y.imag
							T.vmul(d, 2, x + 1, 2, y + 1, 2, n) // y.imag = d.real*x.imag
							T.vma(d + 1, 2, x, 2, y + 1, 2, n)  // y.imag += d.imag*x.real
						}
					}
					complexVMA(A.upper, x + 2, y, &temp, n - 1)
					complexVMA(A.lower, x, y + 2, &temp, n - 1)
				}
			}
		}
	}
	return result
}

/// Matrix–vector multiplication operator.
///
/// Equivalent to `T.multiply(A, x)`.
@inlinable public func *<T: ScalarField>(_ A: TridiagonalMatrix<T>, _ x: [T]) -> [T] { 	T.multiply(A, x) }
