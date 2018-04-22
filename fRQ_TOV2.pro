
CONFIG += c++11

SOURCES += src/main.cpp \
    linalg_library/alglibinternal.cpp \
    linalg_library/alglibmisc.cpp \
    linalg_library/ap.cpp \
    linalg_library/dataanalysis.cpp \
    linalg_library/diffequations.cpp \
    linalg_library/fasttransforms.cpp \
    linalg_library/integration.cpp \
    linalg_library/interpolation.cpp \
    linalg_library/linalg.cpp \
    linalg_library/optimization.cpp \
    linalg_library/solvers.cpp \
    linalg_library/specialfunctions.cpp \
    linalg_library/statistics.cpp \
    src/configure.cpp \
    src/eos_analytical.cpp \
    src/eos_filetoarray.cpp \
    src/eos_functions.cpp \
    src/get_gradient_fr.cpp \
    src/get_gradient_frq.cpp \
    src/get_gradient_gr.cpp \
    src/integrate_star.cpp \
    src/output.cpp \
    src/plotradmassrho.cpp

RESOURCES += qml.qrc

# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =



HEADERS += \
    linalg_library/alglibinternal.h \
    linalg_library/alglibmisc.h \
    linalg_library/ap.h \
    linalg_library/dataanalysis.h \
    linalg_library/diffequations.h \
    linalg_library/fasttransforms.h \
    linalg_library/integration.h \
    linalg_library/interpolation.h \
    linalg_library/linalg.h \
    linalg_library/optimization.h \
    linalg_library/solvers.h \
    linalg_library/specialfunctions.h \
    linalg_library/spline.h \
    linalg_library/statistics.h \
    linalg_library/stdafx.h \
    header/configure.h \
    header/eos_analytical.h \
    header/eos_filetoarray.h \
    header/eos_functions.h \
    header/get_gradient_fr.h \
    header/get_gradient_frq.h \
    header/get_gradient_gr.h \
    header/glob_variables.h \
    header/integrate_star.h \
    header/meta_functions.h \
    header/options.h \
    header/output.h \
    header/plotradmassrho.h

