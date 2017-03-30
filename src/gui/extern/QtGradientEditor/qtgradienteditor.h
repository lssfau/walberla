/****************************************************************************
**
** Copyright (C) 2009 Nokia Corporation and/or its subsidiary(-ies).
** Contact: Qt Software Information (qt-info@nokia.com)
**
** This file is part of the tools applications of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial Usage
** Licensees holding valid Qt Commercial licenses may use this file in
** accordance with the Qt Commercial License Agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Nokia.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** In addition, as a special exception, Nokia gives you certain
** additional rights. These rights are described in the Nokia Qt LGPL
** Exception version 1.0, included in the file LGPL_EXCEPTION.txt in this
** package.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements will be
** met: http://www.gnu.org/copyleft/gpl.html.
**
** If you are unsure which license is appropriate for your use, please
** contact the sales department at qt-sales@nokia.com.
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QTGRADIENTEDITOR_H
#define QTGRADIENTEDITOR_H

#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class QtGradientEditor : public QWidget
{
    Q_OBJECT
    Q_PROPERTY(QGradient gradient READ gradient WRITE setGradient)
    Q_PROPERTY(bool backgroundCheckered READ isBackgroundCheckered WRITE setBackgroundCheckered)
    Q_PROPERTY(bool detailsVisible READ detailsVisible WRITE setDetailsVisible)
    Q_PROPERTY(bool detailsButtonVisible READ isDetailsButtonVisible WRITE setDetailsButtonVisible)
public:
    QtGradientEditor(QWidget *parent = 0);
    ~QtGradientEditor();

    void setGradient(const QGradient &gradient);
    QGradient gradient() const;

    bool isBackgroundCheckered() const;
    void setBackgroundCheckered(bool checkered);

    bool detailsVisible() const;
    void setDetailsVisible(bool visible);

    bool isDetailsButtonVisible() const;
    void setDetailsButtonVisible(bool visible);

    QColor::Spec spec() const;
    void setSpec(QColor::Spec spec);

signals:

    void gradientChanged(const QGradient &gradient);
    void aboutToShowDetails(bool details, int extenstionWidthHint);

private:
    class QtGradientEditorPrivate *d_ptr;
    Q_DECLARE_PRIVATE(QtGradientEditor)
    Q_DISABLE_COPY(QtGradientEditor)
    Q_PRIVATE_SLOT(d_func(), void slotGradientStopsChanged(const QGradientStops &stops))
    Q_PRIVATE_SLOT(d_func(), void slotTypeChanged(int type))
    Q_PRIVATE_SLOT(d_func(), void slotSpreadChanged(int type))
    Q_PRIVATE_SLOT(d_func(), void slotStartLinearXChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotStartLinearYChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotEndLinearXChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotEndLinearYChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotCentralRadialXChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotCentralRadialYChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotFocalRadialXChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotFocalRadialYChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotRadiusRadialChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotCentralConicalXChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotCentralConicalYChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotAngleConicalChanged(double value))
    Q_PRIVATE_SLOT(d_func(), void slotDetailsChanged(bool details))
    Q_PRIVATE_SLOT(d_func(), void startLinearChanged(const QPointF &))
    Q_PRIVATE_SLOT(d_func(), void endLinearChanged(const QPointF &))
    Q_PRIVATE_SLOT(d_func(), void centralRadialChanged(const QPointF &))
    Q_PRIVATE_SLOT(d_func(), void focalRadialChanged(const QPointF &))
    Q_PRIVATE_SLOT(d_func(), void radiusRadialChanged(qreal))
    Q_PRIVATE_SLOT(d_func(), void centralConicalChanged(const QPointF &))
    Q_PRIVATE_SLOT(d_func(), void angleConicalChanged(qreal))
};


#include "qtgradientstopscontroller.h"
#include "ui_qtgradienteditor.h"

#include <QDoubleSpinBox>

class QtGradientEditorPrivate
{
    QtGradientEditor *q_ptr;
    Q_DECLARE_PUBLIC(QtGradientEditor)
public:
    QtGradientEditorPrivate() : m_gradient(QLinearGradient()) {}

    void slotGradientStopsChanged(const QGradientStops &stops);
    void slotTypeChanged(int type);
    void slotSpreadChanged(int spread);
    void slotStartLinearXChanged(double value);
    void slotStartLinearYChanged(double value);
    void slotEndLinearXChanged(double value);
    void slotEndLinearYChanged(double value);
    void slotCentralRadialXChanged(double value);
    void slotCentralRadialYChanged(double value);
    void slotFocalRadialXChanged(double value);
    void slotFocalRadialYChanged(double value);
    void slotRadiusRadialChanged(double value);
    void slotCentralConicalXChanged(double value);
    void slotCentralConicalYChanged(double value);
    void slotAngleConicalChanged(double value);

    void slotDetailsChanged(bool details);

    void startLinearChanged(const QPointF &point);
    void endLinearChanged(const QPointF &point);
    void centralRadialChanged(const QPointF &point);
    void focalRadialChanged(const QPointF &point);
    void radiusRadialChanged(qreal radius);
    void centralConicalChanged(const QPointF &point);
    void angleConicalChanged(qreal angle);

    void setStartLinear(const QPointF &point);
    void setEndLinear(const QPointF &point);
    void setCentralRadial(const QPointF &point);
    void setFocalRadial(const QPointF &point);
    void setRadiusRadial(qreal radius);
    void setCentralConical(const QPointF &point);
    void setAngleConical(qreal angle);

    void setType(QGradient::Type type);
    void showDetails(bool details);

    void setSpinBox(QDoubleSpinBox *spinBox, const char *slot, double max = 1.0, double step = 0.01, int decimals = 3);
    void reset();
    void setLayout(bool details);
    void layoutDetails(bool details);
    bool row4Visible() const;
    bool row5Visible() const;
    int extensionWidthHint() const;

    void setCombos(bool combos);

    QGradient gradient() const;
    void updateGradient(bool emitSignal);

    Ui::QtGradientEditor m_ui;
    QtGradientStopsController *m_gradientStopsController;

    QDoubleSpinBox *startLinearXSpinBox;
    QDoubleSpinBox *startLinearYSpinBox;
    QDoubleSpinBox *endLinearXSpinBox;
    QDoubleSpinBox *endLinearYSpinBox;
    QDoubleSpinBox *centralRadialXSpinBox;
    QDoubleSpinBox *centralRadialYSpinBox;
    QDoubleSpinBox *focalRadialXSpinBox;
    QDoubleSpinBox *focalRadialYSpinBox;
    QDoubleSpinBox *radiusRadialSpinBox;
    QDoubleSpinBox *centralConicalXSpinBox;
    QDoubleSpinBox *centralConicalYSpinBox;
    QDoubleSpinBox *angleConicalSpinBox;

    QButtonGroup *m_typeGroup;
    QButtonGroup *m_spreadGroup;

    QGradient::Type m_type;

    QGridLayout *m_gridLayout;
    QWidget *m_hiddenWidget;
    QGridLayout *m_hiddenLayout;
    bool m_details;
    bool m_detailsButtonVisible;
    bool m_backgroundCheckered;

    QGradient m_gradient;

    bool m_combos;
};


QT_END_NAMESPACE

#endif
