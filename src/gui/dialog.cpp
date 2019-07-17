#include <QDir>
#include <QDateTime>
#include <QFileDialog>
#include <QMessageBox>
#include "dialog.h"
#include "ui_dialog.h"

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog),
    proc_(0)
{
    ui->setupUi(this);

    ui->lineEditWork->setText(QDir::homePath());

    proc_ = new QProcess(this);
    proc_->setProcessChannelMode(QProcess::MergedChannels);
    connect(proc_, SIGNAL(readyReadStandardOutput()), this, SLOT(slot_proc_readyReadStandardOutput()));
    connect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));
}

Dialog::~Dialog()
{
    delete ui;
}

void Dialog::closeEvent(QCloseEvent *e)
{
    if ( isProcessRunning() ) {
        e->ignore();
        return;
    }

    e->accept();
}

void Dialog::on_pushButtonVcf_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose VCF file"));
    if (!fileName.isEmpty())
        ui->lineEditVcf->setText(fileName);
}

void Dialog::on_pushButtonEffect_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose allele effect file"));
    if (!fileName.isEmpty())
        ui->lineEditEffect->setText(fileName);
}

void Dialog::on_pushButtonMap_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose genetic linkage map file"));
    if (!fileName.isEmpty())
        ui->lineEditMap->setText(fileName);
}

void Dialog::on_pushButtonPheno_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose phenotype file"));
    if (!fileName.isEmpty())
        ui->lineEditPheno->setText(fileName);
}

void Dialog::on_pushButtonWeight_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose weights file"));
    if (!fileName.isEmpty())
        ui->lineEditWeight->setText(fileName);
}

void Dialog::on_pushButtonWork_clicked()
{
    QString dir = QFileDialog::getExistingDirectory(this, tr("Set output directory"));
    if (!dir.isEmpty())
        ui->lineEditWork->setText(dir);
}

void Dialog::on_buttonBox_accepted()
{
    QString prog = QDir(QApplication::applicationDirPath()).filePath(QLatin1String("cross"));
#ifdef Q_OS_DARWIN
    if (prog.endsWith(QLatin1String(".app/Contents/MacOS"))) {
        QDir dir(prog);
        dir.cdUp();
        dir.cdUp();
        dir.cdUp();
        prog = exe.absolutePath();
    }
#endif

    QStringList args;

    if ( ! ui->lineEditVcf->text().isEmpty() )
        args << QLatin1String("--vcf") << ui->lineEditVcf->text();

    if ( ! ui->lineEditEffect->text().isEmpty() )
        args << QLatin1String("--effect") << ui->lineEditEffect->text();

    if ( ! ui->lineEditMap->text().isEmpty() )
        args << QLatin1String("--map") << ui->lineEditMap->text();

    if ( ! ui->lineEditPheno->text().isEmpty() )
        args << QLatin1String("--pheno") << ui->lineEditPheno->text();

    if ( ! ui->lineEditWeight->text().isEmpty() )
        args << QLatin1String("--wt") << ui->lineEditWeight->text();

    args << QLatin1String("--type") << QString::number(ui->comboBoxType->currentIndex() + 2);

    if ( ! ui->lineEditSize->text().isEmpty() )
        args << QLatin1String("--size") << ui->lineEditSize->text();

    if ( ! ui->lineEditPct->text().isEmpty() )
        args << QLatin1String("--pct") << ui->lineEditPct->text();

    QString prefix = QLatin1String("cross.out.");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    args << QLatin1String("--out") << QDir(ui->lineEditWork->text()).absoluteFilePath(prefix);

    startProcess(prog, args);
}

void Dialog::on_buttonBox_rejected()
{
    if (ui->buttonBox->button(QDialogButtonBox::Ok)->isEnabled())
        close();

    isProcessRunning();
}

void Dialog::slot_proc_readyReadStandardOutput()
{
    QByteArray v = proc_->readAllStandardOutput();
    ui->plainTextEditLog->moveCursor(QTextCursor::End);
    ui->plainTextEditLog->insertPlainText(QString::fromLocal8Bit(v.data(),v.size()));
    ui->plainTextEditLog->moveCursor(QTextCursor::End);
}

void Dialog::slot_proc_finished(int code, QProcess::ExitStatus status)
{
    ui->progressBar->setRange(0,1);
    ui->progressBar->reset();
    ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);

    if (status != QProcess::NormalExit || code != 0) {
        QMessageBox::critical(this, tr("ERROR"), tr("Process exited unexpectedly: code %1, status %2.").arg(code).arg(status));
        return;
    }
}

bool Dialog::isProcessRunning()
{
    if (proc_->state() != QProcess::NotRunning) {
        if (QMessageBox::question(this, tr("Terminate"), tr("Stop currently running computations?"), QMessageBox::Ok | QMessageBox::Cancel) != QMessageBox::Ok)
            return true;

        disconnect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));
        proc_->close();

        ui->progressBar->setRange(0,1);
        ui->progressBar->reset();
        ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);

        connect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));
    }

    return false;
}

void Dialog::startProcess(const QString &prog, const QStringList &args)
{
    if ( isProcessRunning() )
        return;

    proc_->setWorkingDirectory(ui->lineEditWork->text());

    proc_->start(prog, args);
    if ( ! proc_->waitForStarted() ) {
        QMessageBox::critical(this, tr("ERROR"), tr("Can't start process: %1").arg(prog));
        return;
    }

    ui->progressBar->setRange(0,0);
    ui->progressBar->reset();
    ui->buttonBox->button(QDialogButtonBox::Ok)->setDisabled(true);
}
