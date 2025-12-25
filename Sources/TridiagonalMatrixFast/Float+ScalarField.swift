//
//  Float+ScalarField.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

extension Float: ScalarField, RealScalar {
	
	// MARK: ScalarField (real)
	
	public static var gttrf: GttrfFunction<Float> { sgttrf_ }
	public static var gttrs: GttrsFunction<Float> { sgttrs_ }
	
	public static var gtcon: GtconFunction<Float, Float, Float> {
		{ norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info in
			gtcon_(sgtcon_, norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
		}
	}
	
	public static var axpy: BLASAxpyFunction<Float> {
		{ n, a, x, incx, y, incy in
			axpy_(cblas_saxpy, n, a, x, incx, y, incy)
		}
	}
	
	public static var AXpY: MultiplyAdd<Float> { AXpY_ }
	public static var multiply: MatrixVectorMultiply<Float> { multiply_ }
	
	// MARK: RealScalar
	
	public typealias WType = __CLPK_complex
	
	public static var cgttrf: GttrfFunction<WType> { cgttrf_ }
	public static var cgttrs: GttrsFunction<WType> { cgttrs_ }
	
	public static var cgtcon: GtconFunction<WType, Float, WType> {
		{ norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info in
			gtcon_(cgtcon_, norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
		}
	}
	
	public static var caxpy: BLASAxpyFunction<WType> {
		{ n, a, x, incx, y, incy in
			axpy_(cblas_caxpy, n, a, x, incx, y, incy)
		}
	}
	
	public static var vma: DSPFunction<Float> {
		{ vDSP_vma($0, vDSP_Stride($1), $2, vDSP_Stride($3),
				   $4, vDSP_Stride($5), $4, vDSP_Stride($5),
				   vDSP_Length($6)) }
	}
	
	public static var vmul: DSPFunction<Float> {
		{ vDSP_vmul($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					$4, vDSP_Stride($5), vDSP_Length($6)) }
	}
	
	public static var vsub: DSPFunction<Float> {
		{ vDSP_vsub($0, vDSP_Stride($1), $2, vDSP_Stride($3),
					$4, vDSP_Stride($5), vDSP_Length($6)) }
	}
	
	public static var cAXpY: MultiplyAdd<Complex<Float>> { AXpY_ }
	public static var cMultiply: MatrixVectorMultiply<Complex<Float>> { multiply_ }
}

