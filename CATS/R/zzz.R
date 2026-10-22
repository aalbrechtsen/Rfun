.onLoad=function(libname, pkgname)
{
library.dynam("CATS", pkgname, libname)
}
.onUnload=function(libpath)
{
library.dynam.unload("CATS", libpath)
}