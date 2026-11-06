#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>

template <typename T, std::size_t ROWS, std::size_t COLS>
class Matrix2D {
   private:
	T data[ROWS][COLS];

   public:
	Matrix2D()
	{
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS; ++j) {
				data[i][j] = T{};
			}
		}
	}

	// Инициализация списком списков
	Matrix2D(std::initializer_list<std::initializer_list<T>> init)
	{
		std::size_t i = 0;
		for (auto& rowList : init) {
			std::size_t j = 0;
			for (auto& val : rowList) {
				if (i < ROWS && j < COLS) {
					data[i][j] = val;
				}
				++j;
			}
			++i;
		}
	}

	// Умножение матрицы на скаляр (правый операнд - число)
	Matrix2D<T, ROWS, COLS> operator*(const T& scalar) const
	{
		Matrix2D<T, ROWS, COLS> result;
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS; ++j) {
				result(i, j) = this->operator()(i, j) * scalar;
			}
		}
		return result;
	}

	Matrix2D<T, ROWS, COLS> operator+(const Matrix2D<T, ROWS, COLS>& rhs) const
	{
		Matrix2D<T, ROWS, COLS> result;
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS; ++j) {
				result(i, j) = (*this)(i, j) + rhs(i, j);
			}
		}
		return result;
	}

	// Доступ к элементу (r, c)
	T& operator()(std::size_t r, std::size_t c)
	{
		if (r >= ROWS || c >= COLS) {
			throw std::out_of_range("Index out of range");
		}
		return data[r][c];
	}

	const T& operator()(std::size_t r, std::size_t c) const
	{
		if (r >= ROWS || c >= COLS) {
			throw std::out_of_range("Index out of range");
		}
		return data[r][c];
	}

	// Пример оператора умножения: (A is ROWS×COMMON) * (B is COMMON×COLS) => (RESULT is ROWS×COLS)
	template <std::size_t COMMON, std::size_t COLS2>
	Matrix2D<T, ROWS, COLS2> operator*(const Matrix2D<T, COMMON, COLS2>& rhs) const
	{
		static_assert(COLS == COMMON, "Matrix dimensions mismatch for multiplication");
		Matrix2D<T, ROWS, COLS2> result;
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS2; ++j) {
				T sum = T{};
				for (std::size_t k = 0; k < COMMON; ++k) {
					sum += (*this)(i, k) * rhs(k, j);
				}
				result(i, j) = sum;
			}
		}
		return result;
	}

	// Транспонирование
	Matrix2D<T, COLS, ROWS> transpose() const
	{
		Matrix2D<T, COLS, ROWS> result;
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS; ++j) {
				result(j, i) = data[i][j];
			}
		}
		return result;
	}

	Matrix2D<T, ROWS, COLS> inverse() const
	{
		static_assert(ROWS == COLS, "Inverse can only be computed for square matrices.");
		constexpr std::size_t N = ROWS;

		// Создаём расширенную матрицу [A | I]
		Matrix2D<T, N, 2 * N> augmented;

		// Инициализация расширенной матрицы
		for (std::size_t i = 0; i < N; ++i) {
			for (std::size_t j = 0; j < N; ++j) {
				augmented(i, j) = (*this)(i, j);
			}
			for (std::size_t j = N; j < 2 * N; ++j) {
				augmented(i, j) = (j - N == i) ? static_cast<T>(1) : static_cast<T>(0);
			}
		}

		// Применение метода Гаусса-Жордана
		for (std::size_t i = 0; i < N; ++i) {
			// Поиск главного элемента (максимальный по абсолютной величине в столбце i)
			std::size_t pivotRow = i;
			T maxElem = std::abs(augmented(i, i));
			for (std::size_t k = i + 1; k < N; ++k) {
				if (std::abs(augmented(k, i)) > std::abs(maxElem)) {
					maxElem = std::abs(augmented(k, i));
					pivotRow = k;
				}
			}

			// Проверка на вырожденность матрицы
			if (std::abs(augmented(pivotRow, i)) < static_cast<T>(1e-12)) {
				throw std::runtime_error("Matrix is singular and cannot be inverted.");
			}

			// Поменять текущую строку с pivotRow, если они различны
			if (pivotRow != i) {
				for (std::size_t j = 0; j < 2 * N; ++j) {
					std::swap(augmented(i, j), augmented(pivotRow, j));
				}
			}

			// Нормализовать текущую строку
			T pivot = augmented(i, i);
			for (std::size_t j = 0; j < 2 * N; ++j) {
				augmented(i, j) /= pivot;
			}

			// Обнулить все остальные элементы в текущем столбце
			for (std::size_t k = 0; k < N; ++k) {
				if (k == i)
					continue;
				T factor = augmented(k, i);
				for (std::size_t j = 0; j < 2 * N; ++j) {
					augmented(k, j) -= factor * augmented(i, j);
				}
			}
		}

		// Извлечь обратную матрицу из расширенной матрицы
		Matrix2D<T, ROWS, COLS> inverseMatrix;
		for (std::size_t i = 0; i < N; ++i) {
			for (std::size_t j = 0; j < N; ++j) {
				inverseMatrix(i, j) = augmented(i, j + N);
			}
		}

		return inverseMatrix;
	}

	Matrix2D<T, ROWS, COLS> normalized() const
	{
		T sumSq = T{};
		// Вычисляем сумму квадратов элементов
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS; ++j) {
				sumSq += data[i][j] * data[i][j];
			}
		}
		T norm = std::sqrt(sumSq);
		if (norm == T{}) {
			throw std::runtime_error("Cannot normalize a zero matrix.");
		}

		Matrix2D<T, ROWS, COLS> result;
		// Делим каждый элемент на норму
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS; ++j) {
				result(i, j) = data[i][j] / norm;
			}
		}
		return result;
	}

	// Функция вычисления нормы матрицы (Фробениус)
	T norm() const
	{
		T sumSq = T{};
		for (std::size_t i = 0; i < ROWS; ++i) {
			for (std::size_t j = 0; j < COLS; ++j) {
				sumSq += data[i][j] * data[i][j];
			}
		}
		return std::sqrt(sumSq);
	}
};

#endif  // MATRIX2D_H
