//
//  Double+ScalarField.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

extension Double: ScalarField, RealScalar {
	
	// MARK: ScalarField (real)
	
	public static var gttrf: gttrf<Double> { dgttrf_ }
	public static var gttrs: gttrs<Double> { dgttrs_ }
	
	public static var gtcon: gtcon<Double, Double, Double> {
		{ norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info in
			gtcon_(dgtcon_, norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
		}
	}
	
	public static var axpy: axpy<Double> {
		{ n, a, x, incx, y, incy in
			axpy_(cblas_daxpy, n, a, x, incx, y, incy)
		}
	}
	
	public static var AXpY: MultiplyAdd<Double> { AXpY_ }
	public static var multiply: MatrixVectorMultiply<Double> { multiply_ }
	
	// MARK: RealScalar
	
	public typealias WType = __CLPK_doublecomplex
	
	public static var cgttrf: gttrf<WType> { zgttrf_ }
	public static var cgttrs: gttrs<WType> { zgttrs_ }
	
	public static var cgtcon: gtcon<WType, Double, WType> {
		{ norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info in
			gtcon_(zgtcon_, norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
		}
	}
	
	public static var caxpy: axpy<WType> {
		{ n, a, x, incx, y, incy in
			axpy_(cblas_zaxpy, n, a, x, incx, y, incy)
		}
	}
	
	public static var vma: DSPSignature<Double> {
		{ vDSP_vmaD($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					$4, vDSP_Stride($5), $4, vDSP_Stride($5),
					vDSP_Length($6)) }
	}
	
	public static var vmul: DSPSignature<Double> {
		{ vDSP_vmulD($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					 $4, vDSP_Stride($5), vDSP_Length($6)) }
	}
	
	public static var vsub: DSPSignature<Double> {
		{ vDSP_vsubD($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					 $4, vDSP_Stride($5), vDSP_Length($6)) }
	}
	
	public static var cAXpY: MultiplyAdd<Complex<Double>> { AXpY_ }
	public static var cMultiply: MatrixVectorMultiply<Complex<Double>> { multiply_ }
}
