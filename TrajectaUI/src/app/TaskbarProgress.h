#pragma once

#include <QAbstractNativeEventFilter>
#include <QPointer>

class QWindow;
struct ITaskbarList3;

class TaskbarProgress final : public QAbstractNativeEventFilter {
   public:
	explicit TaskbarProgress(QWindow* window);
	~TaskbarProgress() override;
	void setBusy(bool busy);
	bool nativeEventFilter(const QByteArray& eventType, void* message, qintptr* result) override;

   private:
	void update();
	QPointer<QWindow> m_window;
	ITaskbarList3* m_taskbar{};
	unsigned int m_taskbarCreated{};
	bool m_comInitialized{};
	bool m_busy{};
};
