//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file MainWindow.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BlockSliceView.h"
#include "BlockTreeView.h"
#include "BlockViewText.h"
#include "CMakeDefs.h"
#include "MainWindow.h"

#include "core/DataTypes.h"

#include <QDebug>
#include <QMdiSubWindow>
#include <QMessageBox>
#include <QSettings>


#ifdef WALBERLA_GUI_USE_3D_BLOCKVIEW
#include "BlockView3D.h"
#endif


extern int qInitResources_icons();

namespace walberla {
namespace gui {


MainWindow::MainWindow(timeloop::ITimeloop & timeloop, StructuredBlockForest & blockForest, const GUI & gui)
   : timeloop_(timeloop), blockForest_(blockForest), gui_(gui)
{
   // loading resources that are placed in library
   ::qInitResources_icons();

   ui.setupUi(this);

   ui.cmdContinue->setDefaultAction( ui.actionContinue );
   connect( &timer, SIGNAL(timeout()), this, SLOT(simulationTimerEvent()));


#ifdef WALBERLA_GUI_USE_3D_BLOCKVIEW
   // Block View 3D
   BlockView3D * blockView3D = new BlockView3D();
   sliceChangeListener_ = blockView3D;
   ui.dockBlocks3D->setWidget( blockView3D );
   blockView3D->setup( &blockForest.getBlockForest() );
   connect( this, SIGNAL( blockforestChanged() ), blockView3D, SLOT(onBlockForestChange() ) );
#else
   delete ui.dockBlocks3D;
   sliceChangeListener_ = NULL;
#endif


   // Block View Text
   BlockViewText * blockViewText = new BlockViewText();
   ui.dockBlocksText->setWidget( blockViewText );
   blockViewText->setup ( &blockForest.getBlockForest() );
   connect( this, SIGNAL( blockforestChanged() ), blockViewText, SLOT(onBlockForestChange() ) );


   QSettings settings;
   settings.beginGroup("MainWindowState");
   QByteArray d = settings.value("mainwindow").toByteArray();
   if(d.size() > 0)
     restoreState(d);
   settings.endGroup();

   blockForestModificationStamp_ = blockForest_.getBlockForest().getModificationStamp();
}


void MainWindow::checkForBlockforestChange()
{
   uint_t currentBfStamp = blockForest_.getBlockForest().getModificationStamp();
   if ( blockForestModificationStamp_ != currentBfStamp ) {
      blockForestModificationStamp_ = currentBfStamp;
      emit blockforestChanged();
   }
}


void MainWindow::closeEvent(QCloseEvent * e)
{
   // Destroy BlockSliceView before the BlockView3D
   // -> make sure slice indicators are unregistered before BlockView3D is destroyed
   ui.mdiArea->closeAllSubWindows();

   QSettings settings;
   settings.beginGroup("MainWindowState");
   settings.setValue("mainwindow", saveState());
   settings.endGroup();

   QMainWindow::closeEvent(e);
}


void MainWindow::on_actionNewGridView_triggered()
{
   QWidget *par = new QWidget();

   BlockSliceView * widget = new BlockSliceView( blockForest_, gui_,  ui.propertyStack, sliceChangeListener_ );

   connect(this,   SIGNAL(dataChanged()),
           widget, SLOT( redraw()) );

   QVBoxLayout * boxLayout = new QVBoxLayout(par);
   boxLayout->addWidget( widget );

   QMdiSubWindow *win = ui.mdiArea->addSubWindow( par );
   win->setWindowTitle("Grid View");

   win->resize( QSize(1200,600) );

   par->show();
   widget->show();
}



void MainWindow::on_actionNewTreeView_triggered()
{

   QWidget *par = new QWidget();

   BlockTreeView * widget = new BlockTreeView( gui_ );

   connect(this, SIGNAL(dataChanged()),
           widget, SLOT(onDataChange()) );


   QVBoxLayout * boxLayout = new QVBoxLayout(par);
   boxLayout->addWidget(widget);

   QMdiSubWindow *win = ui.mdiArea->addSubWindow(par);
   win->setWindowTitle("Tree View");

   win->resize(QSize(600,800));

   par->show();
   widget->show();
}

void MainWindow::on_inpTimerInterval_valueChanged(int newValue)
{
    timer.stop();
    timer.start(newValue);
}


void MainWindow::simulationTimerEvent()
{
    simulate(1);
}

void MainWindow::on_actionRun_triggered(bool newState)
{
    if(newState)
    {
        timer.start(ui.inpTimerInterval->value());
        ui.actionPause->setChecked(false);
    }
    else
        ui.actionRun->setChecked(true); //displayed with Pause-action
}

void MainWindow::on_actionPause_triggered(bool newState)
{
    if(newState)
    {
        ui.actionRun->setChecked(false);
        timer.stop();
    }
    else
        ui.actionPause->setChecked(true);
}

void MainWindow::on_actionForward_triggered()
{
    simulate( uint_c( ui.spnTimeStepsForward->value() ) );
}

void MainWindow::on_actionFastForward_triggered()
{
    simulate( uint_c( ui.spnTimeStepsFastForward->value() ) );
}



void MainWindow::simulate(uint_t numTimeSteps)
{
   try
   {
      for(uint_t i=0; i<numTimeSteps; ++i) {
         timeloop_.singleStep();
         ui.lblTimestepsExecuted->setText( QString("%1").arg( timeloop_.getCurrentTimeStep() ) );
         QApplication::instance()->processEvents();
         if ( ui.chkUpdateView->checkState() == Qt::Checked ) {
            emit dataChanged();
            checkForBlockforestChange();
         }
      }
      if ( ui.chkUpdateView->checkState() != Qt::Checked ) {
         emit dataChanged();
         checkForBlockforestChange();
      }
   }
   catch ( std::exception & e )
   {
      QMessageBox::critical( NULL, "Runtime Error", QString::fromStdString( e.what() ) );
      timer.stop();
   }
   checkForBlockforestChange();
}




void MainWindow::breakpoint( const QString & message, const QString & file, int lineNr )
{
   // Check if breakpoint is activated
   if ( ui.chkOnlyTimestep->isChecked() )
   {
      if ( int_c( timeloop_.getCurrentTimeStep() ) != ui.spnOnlyTimestep->value() )
         return;
   }
   if ( ui.chkOnlyLocation->isChecked() )
   {
      if ( ! ui.txtOnlyLocation->text().isEmpty() &&
             file.contains( ui.txtOnlyLocation->text(), Qt::CaseInsensitive ) )
         return;
   }

   ui.actionForward->setEnabled( false );
   ui.actionFastForward->setEnabled( false );
   ui.actionRun->setEnabled( false );
   ui.actionPause->setEnabled( false );

   QPalette pal       = ui.dockBreakpoints->palette();
   pal.setColor(ui.dockBreakpoints->backgroundRole(), Qt::red);
   ui.dockBreakpoints->setPalette(pal);
   ui.dockBreakpoints->raise();


   ui.actionContinue->setEnabled( true );


   ui.lblBreakpointMsg->setText( message );
   ui.lblLocation->setText( QString("%1 : %2").arg(file).arg(lineNr) );


   QEventLoop eventLoop;
   QObject::connect( ui.actionContinue, SIGNAL( triggered() ), &eventLoop, SLOT(quit()) );
   emit dataChanged();
   eventLoop.exec();

   pal       = ui.dockBreakpoints->palette();
   pal.setColor(ui.dockBreakpoints->backgroundRole(), Qt::gray);
   ui.dockBreakpoints->setPalette(pal);


   ui.lblBreakpointMsg->setText( "" );
   ui.lblLocation->setText( "" );
   ui.actionContinue->setEnabled( false );

   ui.actionForward->setEnabled( true );
   ui.actionFastForward->setEnabled( true );
   ui.actionRun->setEnabled( true );
   ui.actionPause->setEnabled( true );
}




} // namespace gui
} // namespace walberla

