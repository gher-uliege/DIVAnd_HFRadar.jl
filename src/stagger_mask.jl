
# Copyright (C) 2009 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.

function stagger_mask(mask_r,op = &)
    sz = [size(mask_r)...]
    sz_u = copy(sz)
    sz_u[1] = sz[1]-1

    sz_v = copy(sz)
    sz_v[2] = sz[2]-1

    sz_psi = copy(sz);
    sz_psi[1] = sz[1]-1
    sz_psi[2] = sz[2]-1


    mask_u = op.(mask_r[1:end-1,:,:],mask_r[2:end,:,:])
    mask_v = op.(mask_r[:,1:end-1,:],mask_r[:,2:end,:])
    mask_psi = op.(mask_r[1:end-1,1:end-1,:],mask_r[1:end-1,2:end,:],mask_r[2:end,1:end-1,:],mask_r[2:end,2:end,:])

    return reshape(mask_u,(sz_u...,)),reshape(mask_v,(sz_v...,)),reshape(mask_psi,(sz_psi...,))
end
