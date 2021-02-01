/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2021 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <stdlib.h>
#include <iostream>

#include "DMRG.h"

using std::cout;
using std::endl;

void CheMPS2::DMRG::PrintLicense(){

   cout << "\n"
           "   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry\n"
           "   Copyright (C) 2013-2021 Sebastian Wouters\n"
           "\n"
           "   This program is free software; you can redistribute it and/or modify\n"
           "   it under the terms of the GNU General Public License as published by\n"
           "   the Free Software Foundation; either version 2 of the License, or\n"
           "   (at your option) any later version.\n"
           "\n"
           "   This program is distributed in the hope that it will be useful,\n"
           "   but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
           "   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
           "   GNU General Public License for more details.\n"
           "\n"
           "   You should have received a copy of the GNU General Public License along\n"
           "   with this program; if not, write to the Free Software Foundation, Inc.,\n"
           "   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.\n" << endl;

}

