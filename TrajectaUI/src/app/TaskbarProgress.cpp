#include "TaskbarProgress.h"

#include <QCoreApplication>
#include <QWindow>
#include <qt_windows.h>
#include <shobjidl.h>

TaskbarProgress::TaskbarProgress(QWindow* window) : m_window(window)
{
	const HRESULT initialized = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
	m_comInitialized = SUCCEEDED(initialized);
	if (SUCCEEDED(CoCreateInstance(CLSID_TaskbarList, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&m_taskbar)))) {
		if (FAILED(m_taskbar->HrInit())) {
			m_taskbar->Release();
			m_taskbar = nullptr;
		}
	}
	m_taskbarCreated = RegisterWindowMessageW(L"TaskbarButtonCreated");
	QCoreApplication::instance()->installNativeEventFilter(this);
}

TaskbarProgress::~TaskbarProgress()
{
	QCoreApplication::instance()->removeNativeEventFilter(this);
	m_busy = false;
	update();
	if (m_taskbar)
		m_taskbar->Release();
	if (m_comInitialized)
		CoUninitialize();
}

void TaskbarProgress::setBusy(bool busy)
{
	if (m_busy == busy)
		return;
	m_busy = busy;
	update();
}

void TaskbarProgress::update()
{
	if (m_taskbar && m_window)
		m_taskbar->SetProgressState(reinterpret_cast<HWND>(m_window->winId()),
		                            m_busy ? TBPF_INDETERMINATE : TBPF_NOPROGRESS);
}

bool TaskbarProgress::nativeEventFilter(const QByteArray&, void* message, qintptr*)
{
	const auto* event = static_cast<const MSG*>(message);
	if (m_taskbarCreated && event && event->message == m_taskbarCreated && m_window &&
	    event->hwnd == reinterpret_cast<HWND>(m_window->winId()))
		update();
	return false;
}
