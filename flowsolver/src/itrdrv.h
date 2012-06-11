#ifndef _ITRDRV_H_
#define _ITRDRV_H_

extern "C" void itrdrv_init();
extern "C" void itrdrv_iter_init();
extern "C" void itrdrv_iter_step();
extern "C" void itrdrv_iter_advance();
extern "C" void itrdrv_iter_finalize();
extern "C" void itrdrv_finalize();

#endif
