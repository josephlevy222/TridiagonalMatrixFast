//
//  Complex+ScalarField.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

extension Complex: ScalarField where RealType: RealScalar {
	
	public typealias CType = RealType.WType
	
	public static var gttrf: gttrf<Complex<RealType>> {
		{ n, dl, d, du, du2, ipiv, info in
			reboundBands(
				n: Int(n!.pointee),
				dl: dl!,
				d: d!,
				du: du!,
				du2: du2!,
				to: CType.self
			) { dlC, dC, duC, du2C in
				RealType.cgttrf(n, dlC, dC, duC, du2C, ipiv, info)
			}
		}
	}
	
	public static var gttrs: gttrs<Complex<RealType>> {
		{ trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info in
			guard let b else { return -1 }
			let countN = Int(n!.pointee)
			let countB = max(Int(nrhs!.pointee) * countN, 1)
			return b.withMemoryRebound(to: CType.self, capacity: countB) { bC in
				reboundBands(
					n: countN,
					dl: dl!,
					d: d!,
					du: du!,
					du2: du2!,
					to: CType.self
				) { dlC, dC, duC, du2C in
					RealType.cgttrs(
						trans, n, nrhs,
						dlC, dC, duC, du2C,
						ipiv, bC, ldb, info
					)
				}
			}
		}
	}
	
	public static var gtcon: gtcon<Complex<RealType>, RealType, CType> {
		{ norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info in
			reboundBands(
				n: Int(n!.pointee),
				dl: dl!,
				d: d!,
				du: du!,
				du2: du2!,
				to: CType.self
			) { dlC, dC, duC, du2C in
				RealType.cgtcon(
					norm, n,
					dlC, dC, duC, du2C,
					ipiv, anorm, rcond, work, iwork, info
				)
			}
		}
	}
	
	public static var axpy: axpy<Complex<RealType>> {
		{ n, a, x, incx, y, incy in
			withUnsafeBytes(of: a) { aRaw in
				let aC = aRaw.bindMemory(to: CType.self).baseAddress!
				x.withMemoryRebound(to: CType.self, capacity: Int(n)) { xC in
					y.withMemoryRebound(to: CType.self, capacity: Int(n)) { yC in
						RealType.caxpy(n, aC.pointee, xC, incx, yC, incy)
					}
				}
			}
		}
	}
	
	public static var AXpY: MultiplyAdd<Complex<RealType>> { RealType.cAXpY }
	public static var multiply: MatrixVectorMultiply<Complex<RealType>> { RealType.cMultiply }
	
	@inline(__always)
	static func reboundBands<T, CType, R>(
		n: Int,
		dl: CMutablePtr<T>,
		d: CMutablePtr<T>,
		du: CMutablePtr<T>,
		du2: CMutablePtr<T>,
		to: CType.Type = CType.self,
		body: (CMutablePtr<CType>, CMutablePtr<CType>, CMutablePtr<CType>, CMutablePtr<CType>) -> R
	) -> R {
		precondition(n > 0, "n must be > 0")
		return dl.withMemoryRebound(to: CType.self, capacity: max(1, n - 1)) { dlC in
			d.withMemoryRebound(to: CType.self, capacity: n) { dC in
				du.withMemoryRebound(to: CType.self, capacity: max(1, n - 1)) { duC in
					du2.withMemoryRebound(to: CType.self, capacity: max(0, n - 2)) { du2C in
						body(dlC, dC, duC, du2C)
					}
				}
			}
		}
	}
}
