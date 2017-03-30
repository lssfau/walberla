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
//! \file MainWindow.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef Q_MOC_RUN

#include "Gui.h"
#include "ISliceChangeListener.h"
#include "ui_MainWindow.h"

#include "blockforest/StructuredBlockForest.h"

#include "field/GhostLayerField.h"
#include "timeloop/ITimeloop.h"
#endif

#include <QListWidget>
#include <QMainWindow>
#include <QMdiSubWindow>
#include <QTimer>


namespace walberla {
namespace gui {


   //*******************************************************************************************************************
   /*!
   * The main window of the GUI
   *
   * \ingroup gui
   *
   * To setup the window the MainWindow.ui file is used.
   */
   //*******************************************************************************************************************
   class MainWindow : public QMainWindow
   {
   Q_OBJECT

   public:
      MainWindow(timeloop::ITimeloop & timeloop, StructuredBlockForest & blockForest, const GUI & gui );

      virtual ~MainWindow() {}

   signals:
      /// Emitted after simulation progressed and updated data are available
      void dataChanged();
      void blockforestChanged();

   public slots:
      void simulate( uint_t numTimeSteps );

      void breakpoint( const QString & message, const QString & file, int lineNr );

   private slots:

      void simulationTimerEvent();

      // The following slots have "magic" names so that they are automatically connected
      // by the "ui.setupUi(this)" call
      void on_actionNewGridView_triggered();
      void on_actionNewTreeView_triggered();
      void on_actionRun_triggered(bool);
      void on_actionPause_triggered(bool);
      void on_actionForward_triggered();
      void on_actionFastForward_triggered();
      void on_inpTimerInterval_valueChanged(int);

   protected:
      virtual void closeEvent( QCloseEvent * e);

   private:
      void checkForBlockforestChange();

      Ui::MainWindow ui;
      ISliceChangeListener * sliceChangeListener_;

      QTimer timer;

      timeloop::ITimeloop   & timeloop_;
      StructuredBlockForest & blockForest_;
      uint_t                  blockForestModificationStamp_;
      const GUI & gui_;
   };

} // namespace gui
} // namespace walberla


