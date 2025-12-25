//
//  TridiagonalLUMatrix.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

/// The LU factorization of a tridiagonal matrix.
///
/// Stores:
/// - L (unit lower bidiagonal)
/// - U (upper bidiagonal + second superdiagonal)
/// - Pivot indices
/// - Condition estimate
/// - Determinant
@frozen
public struct TridiagonalLUMatrix<T: ScalarField> {
	
	public var lower: ContiguousArray<T>
	public var diagonal: ContiguousArray<T>
	public var upper: ContiguousArray<T>
	public var upper2: ContiguousArray<T>
	public var ipiv: ContiguousArray<CInt>
	
	public var rcond: T.Magnitude
	public var anorm: T.Magnitude
	public var isSingular: Bool
	public var determinant: T
	
	/// Number of rows/columns.
	public var count: Int { diagonal.count }
	
	/// Workspace used for condition estimation.
	private let workspace: TridiagonalWorkspace<T>
	
	/// Approximate condition number `cond(A) ≈ 1 / rcond`.
	public var approximateConditionNumber: T.Magnitude {
		(!isSingular && rcond > 0) ? 1 / rcond : T.Magnitude.infinity
	}
	
	@inline(__always)
	mutating func withMutableLUBufferPointers<R>(
		_ body: (
			_ dl: CMutablePtr<T>,
			_ d: CMutablePtr<T>,
			_ du: CMutablePtr<T>,
			_ du2: CMutablePtr<T>,
			_ ipiv: CMutablePtr<CInt>
		) -> R
	) -> R {
		lower.withUnsafeMutableBufferPointer { dlBuf in
			diagonal.withUnsafeMutableBufferPointer { dBuf in
				upper.withUnsafeMutableBufferPointer { duBuf in
					upper2.withUnsafeMutableBufferPointer { du2Buf in
						ipiv.withUnsafeMutableBufferPointer { ipivBuf in
							guard let dlBase = dlBuf.baseAddress,
								  let dBase = dBuf.baseAddress,
								  let duBase = duBuf.baseAddress,
								  let du2Base = du2Buf.baseAddress,
								  let ipivBase = ipivBuf.baseAddress
							else {
								preconditionFailure("Unexpected nil base address")
							}
							return body(dlBase, dBase, duBase, du2Base, ipivBase)
						}
					}
				}
			}
		}
	}
	
	/// Creates an LU factorization of the given tridiagonal matrix.
	///
	/// - Parameter A: A tridiagonal matrix.
	public init(_ A: TridiagonalMatrix<T>) {
		self.init(dl: A.lower, d: A.diagonal, du: A.upper, anorm: A.oneNorm())
	}
	
	@inline(__always)
	init(dl: [T], d: [T], du: [T], anorm: T.Magnitude) {
		precondition(d.count > 0, "Matrix must be non-empty")
		precondition(dl.count == d.count - 1 && du.count == d.count - 1, "Diagonal sizes inconsistent")
		
		self.lower = ContiguousArray(dl)
		self.diagonal = ContiguousArray(d)
		self.upper = ContiguousArray(du)
		self.upper2 = ContiguousArray(repeating: .zero, count: max(0, d.count - 2))
		self.ipiv = ContiguousArray(repeating: 0, count: d.count)
		self.anorm = anorm
		self.rcond = 0
		self.isSingular = true
		self.determinant = 0
		self.workspace = TridiagonalWorkspace(capacity: d.count)
		
		var n_ = CInt(d.count)
		var info = CInt(0)
		
		// Factor
		withMutableLUBufferPointers { dlPtr, dPtr, duPtr, du2Ptr, ipivPtr in
			_ = T.gttrf(&n_, dlPtr, dPtr, duPtr, du2Ptr, ipivPtr, &info)
		}
		if info != 0 { return }
		isSingular = false
		
		// Determinant: product of U diagonal times permutation parity
		let detU = diagonal.reduce(T.one, *)
		let n = count
		var seen = Array(repeating: false, count: n)
		var cycles = 0
		for i in 0..<n where !seen[i] {
			var j = i
			cycles += 1
			while !seen[j] {
				seen[j] = true
				let pivot = Int(ipiv[j]) - 1
				j = pivot
			}
		}
		let parityIsNegative = ((n - cycles) % 2) != 0
		determinant = parityIsNegative ? -detU : detU
		
		// Condition estimate (gtcon)
		var anorm_ = anorm
		var rcond_ = T.Magnitude(0)
		var normChar: Int8 = Int8(UnicodeScalar("O").value)
		let workPtr = workspace.workBuffer(for: count)
		let iworkPtr = workspace.iworkBuffer(for: count)
		withMutableLUBufferPointers { dlPtr, dPtr, duPtr, du2Ptr, ipivPtr in
			_ = T.gtcon(
				&normChar, &n_,
				dlPtr, dPtr, duPtr, du2Ptr,
				ipivPtr, &anorm_, &rcond_,
				workPtr, iworkPtr, &info
			)
		}
		self.rcond = rcond_
	}
	
	/// Solves `A·x = b` for a single right-hand side.
	///
	/// - Parameters:
	///   - b: Right-hand side (overwritten with the solution).
	///   - transpose: If `true`, solves `Aᵗ·x = b` (real) or `Aᴴ·x = b` (complex).
	/// - Returns: The solution vector.
	@discardableResult
	public mutating func solve(_ b: inout [T], transpose: Bool = false) -> [T] {
		precondition(b.count == self.count)
		guard !isSingular else { return Array(repeating: .zero, count: count) }
		
		var n = CInt(count)
		var nrhs = CInt(1)
		var ldb = n
		var info = CInt(0)
		var trans: Int8 = transpose
		? (T.self is (any RealScalar.Type)
		   ? Int8(UnicodeScalar("T").value)
		   : Int8(UnicodeScalar("C").value))
		: Int8(UnicodeScalar("N").value)
		
		b.withUnsafeMutableBufferPointer { buffer in
			guard let bPtr = buffer.baseAddress else {
				preconditionFailure("Array base address is nil")
			}
			withMutableLUBufferPointers { dlPtr, dPtr, duPtr, du2Ptr, ipivPtr in
				_ = T.gttrs(
					&trans, &n, &nrhs,
					dlPtr, dPtr, duPtr, du2Ptr,
					ipivPtr, bPtr, &ldb, &info
				)
			}
		}
		return b
	}
	
	/// Solves `A·X = B` where `B` is column-major.
	///
	/// - Parameters:
	///   - bColumnMajor: Right-hand sides, column-major (`n * nrhs` elements).
	///   - nrhs: Number of right-hand sides.
	///   - transpose: If `true`, solves `Aᵗ·X = B` (real) or `Aᴴ·X = B` (complex).
	/// - Returns: The solution in the same layout.
	@discardableResult
	public mutating func solve(
		_ bColumnMajor: inout [T],
		nrhs: Int,
		transpose: Bool = false
	) -> [T] {
		precondition(
			bColumnMajor.count == self.count * nrhs,
			"B must be column-major n x nrhs"
		)
		guard !isSingular else {
			return Array(repeating: .zero, count: bColumnMajor.count)
		}
		
		var n = CInt(count)
		var nrhs_c = CInt(nrhs)
		var ldb = n
		var info = CInt(0)
		var trans: Int8 = transpose
		? (T.self is (any RealScalar.Type)
		   ? Int8(UnicodeScalar("T").value)
		   : Int8(UnicodeScalar("C").value))
		: Int8(UnicodeScalar("N").value)
		
		bColumnMajor.withUnsafeMutableBufferPointer { buf in
			guard let bPtr = buf.baseAddress else {
				preconditionFailure("B base address is nil")
			}
			withMutableLUBufferPointers { dlPtr, dPtr, duPtr, du2Ptr, ipivPtr in
				_ = T.gttrs(
					&trans, &n, &nrhs_c,
					dlPtr, dPtr, duPtr, du2Ptr,
					ipivPtr, bPtr, &ldb, &info
				)
			}
		}
		return bColumnMajor
	}
	
	/// Solves `A·X = B` where `B` is an array of column vectors.
	///
	/// - Parameters:
	///   - columns: Right-hand side vectors, length `n` each.
	///   - transpose: If `true`, solves `Aᵗ·X = B` (real) or `Aᴴ·X = B` (complex).
	/// - Returns: The solution columns in the same shape.
	@discardableResult
	public mutating func solve(
		columns: [[T]],
		transpose: Bool = false
	) -> [[T]] {
		let nrhs = columns.count
		guard nrhs > 0 else { return [] }
		let n = count
		precondition(columns.allSatisfy { $0.count == n }, "All columns must have length n")
		
		// Pack into column-major contiguous buffer
		var bFlat = [T](unsafeUninitializedCapacity: n * nrhs) { buffer, initializedCount in
			var offset = 0
			for col in columns {
				col.withUnsafeBufferPointer { colPtr in
					buffer.baseAddress!.advanced(by: offset)
						.initialize(from: colPtr.baseAddress!, count: n)
				}
				offset += n
			}
			initializedCount = n * nrhs
		}
		
		let solvedFlat = solve(&bFlat, nrhs: nrhs, transpose: transpose)
		
		// Unpack
		var out = [[T]]()
		out.reserveCapacity(nrhs)
		solvedFlat.withUnsafeBufferPointer { flatPtr in
			for j in 0..<nrhs {
				let start = j * n
				let col = Array(unsafeUninitializedCapacity: n) { buffer, count in
					buffer.baseAddress!.initialize(
						from: flatPtr.baseAddress! + start,
						count: n
					)
					count = n
				}
				out.append(col)
			}
		}
		return out
	}
}

