/****************************************************************************
**
** Copyright (C) 2011 Andrey Kartashov .
** All rights reserved.
** Contact: Andrey Kartashov (porter@porter.st)
**
** This file is part of the averagedensity module of the genome-tools.
**
** GNU Lesser General Public License Usage
** This file may be used under the terms of the GNU Lesser General Public
** License version 2.1 as published by the Free Software Foundation and
** appearing in the file LICENSE.LGPL included in the packaging of this
** file. Please review the following information to ensure the GNU Lesser
** General Public License version 2.1 requirements will be met:
** http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** Other Usage
** Alternatively, this file may be used in accordance with the terms and
** conditions contained in a signed written agreement between you and Andrey Kartashov.
**
****************************************************************************/

#ifndef _ATDHEATMAP_
#define _ATDHEATMAP_

#include <config.hpp>
#include <atdpbasics.hpp>
#include <Math.hpp>

#ifndef FSTM
#define FSTM ATDP
#endif

class ATDP: public QObject
{
        Q_OBJECT
    private:

        QMap<QString,EXPERIMENT_INFO> experiment_info;
        int avd_window;
        int avd_whole_region;
        int avd_heat_window;
        int avd_bodysize;
        QString twicechr;
        QString ignorechr;

        void getRecordInfo(void);
        void print(ostream& output_stream);
        bool export_to_file(const string & output_filename);
    public slots:
        void start(void);

    signals:
        void finished(void);

    public:

        ATDP(QObject* parent=0);
        ~ATDP();
};

#endif
