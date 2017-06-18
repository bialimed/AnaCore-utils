#
# Copyright (C) 2017 IUCT-O
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class CigarlineGraph (Component):

    def define_parameters(self, aln_files):
        self.add_input_file_list( "aln_files", "Which reads files should be used. ***************************", default=aln_files, required=True )
        self.add_output_file_list( "cigar_R1", "******************", pattern='{basename_woext}_cigarR1.tsv', items=self.aln_files )
        self.add_output_file_list( "cigar_R2", "******************", pattern='{basename_woext}_cigarR2.tsv', items=self.aln_files )
        self.add_output_file_list( "stderr", "The cigarlinegraph stderr file", pattern='{basename_woext}.stderr', items=self.aln_files )

    def process(self):
        cigarlinegraph = ShellFunction( self.get_exec_path("samtools") + ' view -F256 $1 | ' + self.get_exec_path("cigarlineGraph") + " --readsplit -i - -t $2 $3 2> $4", cmd_format='{EXE} {IN} {OUT}' )
        cigarlinegraph = MultiMap( cigarlinegraph, inputs=[self.aln_files], outputs=[self.cigar_R1, self.cigar_R2, self.stderr] )
