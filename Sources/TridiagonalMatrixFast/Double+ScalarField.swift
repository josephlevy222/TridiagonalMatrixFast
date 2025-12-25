//
//  Double+ScalarField.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

extension Double: ScalarField, RealScalar {
	
	// MARK: ScalarField (real)
	
	public static var gttrf: GttrfFunction<Double> { dgttrf_ }
	public static var gttrs: GttrsFunction<Double> { dgttrs_ }
	
	public static var gtcon: GtconFunction<Double, Double, Double> {
		{ norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info in
			gtcon_(dgtcon_, norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
		}
	}
	
	public static var axpy: BLASAxpyFunction<Double> {
		{ n, a, x, incx, y, incy in
			axpy_(cblas_daxpy, n, a, x, incx, y, incy)
		}
	}
	
	public static var AXpY: MultiplyAdd<Double> { AXpY_ }
	public static var multiply: MatrixVectorMultiply<Double> { multiply_ }
	
	// MARK: RealScalar
	
	public typealias WType = __CLPK_doublecomplex
	
	public static var cgttrf: GttrfFunction<WType> { zgttrf_ }
	public static var cgttrs: GttrsFunction<WType> { zgttrs_ }
	
	public static var cgtcon: GtconFunction<WType, Double, WType> {
		{ norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info in
			gtcon_(zgtcon_, norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
		}
	}
	
	public static var caxpy: BLASAxpyFunction<WType> {
		{ n, a, x, incx, y, incy in
			axpy_(cblas_zaxpy, n, a, x, incx, y, incy)
		}
	}
	
	public static var vma: DSPFunction<Double> {
		{ vDSP_vmaD($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					$4, vDSP_Stride($5), $4, vDSP_Stride($5),
					vDSP_Length($6)) }
	}
	
	public static var vmul: DSPFunction<Double> {
		{ vDSP_vmulD($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					 $4, vDSP_Stride($5), vDSP_Length($6)) }
	}
	
	public static var vsub: DSPFunction<Double> {
		{ vDSP_vsubD($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					 $4, vDSP_Stride($5), vDSP_Length($6)) }
	}
	
	public static var cAXpY: MultiplyAdd<Complex<Double>> { AXpY_ }
	public static var cMultiply: MatrixVectorMultiply<Complex<Double>> { multiply_ }
}
