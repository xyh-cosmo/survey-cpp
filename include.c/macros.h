#ifndef _MACROS_H_
#define _MACROS_H_

#define print_Error_ID(rank) \
{if(rank==0) printf("\n==> File: %s\n==> Func: %s\n==> Line: %d\n",__FILE__,__FUNCTION__,__LINE__);}

#define print_Error_Msg(msg,rank) \
{print_Error_ID(rank);\
if(rank==0) printf("Info: %s\n\n",msg);}

#endif