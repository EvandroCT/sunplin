/*************************************************************************
	
	Copyright (C) 2016	Evandro Taquary, Mateus Freitas
	
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*************************************************************************/

/* GPU modular copy */
__global__ void modcpy(void *destination, void *source, size_t destination_size, size_t source_size){

	int idx = blockIdx.x * blockDim.x + threadIdx.x; 
	int pos;
 
	int ds = destination_size/sizeof(int4), ss = source_size/sizeof(int4);
	for(int i = idx; i < ds; i += gridDim.x * blockDim.x){
		pos = i % ss;
		reinterpret_cast<int4*>(destination)[i] = reinterpret_cast<int4*>(source)[pos];  
	}
}
