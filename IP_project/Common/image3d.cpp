#include "image3d.h"
#include "memory.h"
#include <cassert>

using mc::image3d;

/*
* Explicit Template Instantiation.
*/
template class __declspec(dllexport) image3d<char>;
template class __declspec(dllexport) image3d<unsigned char>;
template class __declspec(dllexport) image3d<short>;
template class __declspec(dllexport) image3d<unsigned short>;
template class __declspec(dllexport) image3d<int>;
template class __declspec(dllexport) image3d<unsigned int>;
template class __declspec(dllexport) image3d<float>;

template <typename TYPE>
image3d<TYPE>::image3d(const int w, const int h, const int d) noexcept :
m_ppData(nullptr),
m_iWidth(w),
m_iHeight(h),
m_iDepth(d),
m_fPixelSpacing(1.0f),
m_fSliceSpacing(1.0f)
{
	assert(w != 0 && h != 0 && d != 0);
	m_ppData = SAFE_ALLOC_VOLUME(TYPE, d, w*h);
}

template <typename TYPE>
image3d<TYPE>::~image3d() noexcept
{
	SAFE_DELETE_VOLUME(m_ppData, m_iDepth);
}

template <typename TYPE>
image3d<TYPE>::image3d(const image3d& rhs) noexcept :
m_ppData(nullptr),
m_iWidth(rhs.m_iWidth),
m_iHeight(rhs.m_iHeight),
m_iDepth(rhs.m_iDepth),
m_fPixelSpacing(rhs.m_fPixelSpacing),
m_fSliceSpacing(rhs.m_fSliceSpacing)
{
	m_ppData = SAFE_ALLOC_VOLUME(TYPE, m_iDepth, m_iWidth * m_iHeight);
	MEMCPY_VOLUME(TYPE, m_ppData, rhs.data(), m_iDepth, m_iWidth*m_iHeight);
}

template <typename TYPE>
image3d<TYPE>& image3d<TYPE>::operator=(const image3d& rhs) noexcept
{
	SAFE_DELETE_VOLUME(m_ppData, m_iDepth);

	m_iWidth = rhs.m_iWidth;
	m_iHeight = rhs.m_iHeight;
	m_iDepth = rhs.m_iDepth;
	m_fPixelSpacing = rhs.m_fPixelSpacing;
	m_fSliceSpacing = rhs.m_fSliceSpacing;

	m_ppData = SAFE_ALLOC_VOLUME(TYPE, m_iDepth, m_iWidth * m_iHeight);
	MEMCPY_VOLUME(TYPE, m_ppData, rhs.m_ppData, m_iDepth, m_iWidth * m_iHeight);

	return *this;
}

template <typename TYPE>
image3d<TYPE>::image3d(image3d&& rhs) noexcept :
m_ppData(rhs.m_ppData),
m_iWidth(rhs.m_iWidth),
m_iHeight(rhs.m_iHeight),
m_iDepth(rhs.m_iDepth),
m_fPixelSpacing(rhs.m_fPixelSpacing),
m_fSliceSpacing(rhs.m_fSliceSpacing)
{
	rhs.m_ppData = nullptr;
}

template <typename TYPE>
image3d<TYPE>& image3d<TYPE>::operator=(image3d&& rhs) noexcept
{
	SAFE_DELETE_VOLUME(m_ppData, m_iDepth);

	m_iWidth = rhs.width();
	m_iHeight = rhs.height();
	m_iDepth = rhs.depth();
	m_fPixelSpacing = rhs.pixelSpacing();
	m_fSliceSpacing = rhs.sliceSpacing();

	m_ppData = rhs.m_ppData;
	rhs.m_ppData = nullptr;

	return *this;
}

template <typename TYPE>
void image3d<TYPE>::zeroImage() noexcept
{
	ZERO_VOLUME(TYPE, m_ppData, m_iDepth, this->sizeSlice());
}