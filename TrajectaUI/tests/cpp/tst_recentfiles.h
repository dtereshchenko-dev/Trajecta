#ifndef TST_RECENTFILES_H
#define TST_RECENTFILES_H

#include <QTest>
#include <QSignalSpy>
#include <QSettings>
#include <QTemporaryDir>
#include <QDir>
#include <QDateTime>
#include "io/RecentFiles.h"

class TestRecentFiles : public QObject {
	Q_OBJECT

   private:
	QString createTempFile(const QString& name)
	{
		QString path = m_files.filePath(name);
		QFile f(path);
		if (!f.open(QIODevice::WriteOnly) || f.write("test") != 4) {
			QTest::qFail("Could not create recent file fixture", __FILE__, __LINE__);
			return {};
		}
		f.close();
		m_tempFiles.append(path);
		return path;
	}

	void cleanupTempFiles()
	{
		for (const auto& f : std::as_const(m_tempFiles))
			QFile::remove(f);
		m_tempFiles.clear();
	}

	QStringList m_tempFiles;
	QTemporaryDir m_files;

   private slots:

	void init()
	{
		QVERIFY(m_files.isValid());
		QSettings settings;
		settings.clear();
	}

	void cleanup()
	{
		QSettings settings;
		settings.clear();
		cleanupTempFiles();
	}

	void addFile_increasesCount()
	{
		RecentFiles rf;
		QCOMPARE(rf.count(), 0);
		QString path = createTempFile("test1.txt");
		rf.addFile(QUrl::fromLocalFile(path));
		QCOMPARE(rf.count(), 1);
	}

	void addFile_emitsSignal()
	{
		RecentFiles rf;
		QSignalSpy spy(&rf, &RecentFiles::filesChanged);
		QString path = createTempFile("test1.txt");
		rf.addFile(QUrl::fromLocalFile(path));
		QCOMPARE(spy.count(), 1);
	}

	void addFile_deduplicates()
	{
		RecentFiles rf;
		QString path = createTempFile("test1.txt");
		rf.addFile(QUrl::fromLocalFile(path));
		rf.addFile(QUrl::fromLocalFile(path));
		QCOMPARE(rf.count(), 1);
	}

	void addFile_movesToTop()
	{
		RecentFiles rf;
		QString p1 = createTempFile("first.txt");
		QString p2 = createTempFile("second.txt");
		rf.addFile(QUrl::fromLocalFile(p1));
		rf.addFile(QUrl::fromLocalFile(p2));
		rf.addFile(QUrl::fromLocalFile(p1));
		QCOMPARE(rf.filePath(0), p1);
		QCOMPARE(rf.filePath(1), p2);
	}

	void addFile_maxRecent()
	{
		RecentFiles rf;
		QStringList paths;
		for (int i = 0; i < 7; ++i) {
			paths << createTempFile(QString("file%1.txt").arg(i));
			rf.addFile(QUrl::fromLocalFile(paths.last()));
		}
		QCOMPARE(rf.count(), 6);
		QVERIFY(!rf.files().contains(paths[0]));
		QCOMPARE(rf.filePath(0), paths[6]);
	}

	void addFile_ignoresEmptyUrl()
	{
		RecentFiles rf;
		rf.addFile(QUrl());
		QCOMPARE(rf.count(), 0);
	}

	void removeFile_decreasesCount()
	{
		RecentFiles rf;
		QString path = createTempFile("test1.txt");
		rf.addFile(QUrl::fromLocalFile(path));
		rf.removeFile(path);
		QCOMPARE(rf.count(), 0);
	}

	void removeFile_emitsSignal()
	{
		RecentFiles rf;
		QString path = createTempFile("test1.txt");
		rf.addFile(QUrl::fromLocalFile(path));
		QSignalSpy spy(&rf, &RecentFiles::filesChanged);
		rf.removeFile(path);
		QCOMPARE(spy.count(), 1);
	}

	void removeFile_nonExistent_noSignal()
	{
		RecentFiles rf;
		QSignalSpy spy(&rf, &RecentFiles::filesChanged);
		rf.removeFile("/nonexistent/path.txt");
		QCOMPARE(spy.count(), 0);
	}

	void clear_emptiesList()
	{
		RecentFiles rf;
		rf.addFile(QUrl::fromLocalFile(createTempFile("a.txt")));
		rf.addFile(QUrl::fromLocalFile(createTempFile("b.txt")));
		rf.clear();
		QCOMPARE(rf.count(), 0);
	}

	void fileName_returnsBaseName()
	{
		RecentFiles rf;
		QString path = createTempFile("myfile.mpc");
		rf.addFile(QUrl::fromLocalFile(path));
		QCOMPARE(rf.fileName(0), "myfile.mpc");
	}

	void filePath_returnsFullPath()
	{
		RecentFiles rf;
		QString path = createTempFile("myfile.mpc");
		rf.addFile(QUrl::fromLocalFile(path));
		QCOMPARE(rf.filePath(0), path);
	}

	void displayPath_containsSeparators()
	{
		RecentFiles rf;
		QString path = createTempFile("myfile.mpc");
		rf.addFile(QUrl::fromLocalFile(path));
		QString dp = rf.displayPath(0);
		QVERIFY(!dp.isEmpty());
#ifdef Q_OS_WIN
		QVERIFY(dp.contains('\\'));
#endif
	}

	void getters_invalidIndex_returnEmpty()
	{
		RecentFiles rf;
		QCOMPARE(rf.fileName(-1), QString());
		QCOMPARE(rf.fileName(0), QString());
		QCOMPARE(rf.filePath(99), QString());
		QCOMPARE(rf.displayPath(-5), QString());
	}

	void persistence_savesAndLoads()
	{
		QString path = createTempFile("persist.txt");
		{
			RecentFiles rf;
			rf.addFile(QUrl::fromLocalFile(path));
		}
		{
			RecentFiles rf2;
			QCOMPARE(rf2.count(), 1);
			QCOMPARE(rf2.filePath(0), path);
		}
	}
};

#endif  // TST_RECENTFILES_H
