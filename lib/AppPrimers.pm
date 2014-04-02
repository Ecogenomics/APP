#!/usr/bin/perl
###############################################################################
#
#    AppPrimers.pl
#    
#    Listing of primers used in pyrodb
#
#    Copyright (C) 2011 Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
package AppPrimers;
require Exporter;

our @ISA = qw(Exporter);
our @EXPORT=qw(
      %APP_prim_len_hash
    );


our %APP_prim_len_hash = ();
$APP_prim_len_hash{'pyroL803Fmix'} = 688;
$APP_prim_len_hash{'pyroL926F'} = 565;
$APP_prim_len_hash{'1114'} = 500;
$APP_prim_len_hash{'pyroLSSU926F'} = 565;
$APP_prim_len_hash{'pyroLSSU803F'} = 688;
$APP_prim_len_hash{'pyroLSSU1114F-1'} = 379;
$APP_prim_len_hash{'pyroLSSU1114F-3'} = 379;
$APP_prim_len_hash{'pyroL_ITS1-F'} = 425;
$APP_prim_len_hash{'pyroLfimH72F'} = 640;
$APP_prim_len_hash{'pyroLfimH563R'} = 640;
$APP_prim_len_hash{'ppk_pyro_968'} = 470;
$APP_prim_len_hash{'iTAG926F'} = 580;

1;
