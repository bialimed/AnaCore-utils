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

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from jflow.component import Component
from jflow.abstraction import MultiMap

from weaver.function import ShellFunction


class FilterVCFAnnotOnRNA (Component):

    def define_parameters(self, in_RNA, in_variants, without_version=False):
        # Parameters
        self.add_parameter( "without_version", "With True the version number is removed from the RNA_id.", default=without_version, type=bool )

        # Input files
        self.add_input_file( "in_RNA", "The path to the file describing the RNA kept for each gene (format: TSV). Except the lines starting with a sharp each line has the following format: <GENE>\t<RNA_ID>.", default=in_RNA, required=True )
        self.add_input_file_list( "in_variants", "The path to the variants file (format: VCF).", default=in_variants, required=True )

        # Output files
        self.add_output_file_list( "out_variants", "The path to the outputted variants file (format: VCF).", pattern='{basename_woext}.vcf', items=self.in_variants )
        self.add_output_file_list( "stderr", "The path to the stderr file (format: txt).", pattern='{basename_woext}.stderr', items=self.in_variants )

    def process(self):
        cmd = self.get_exec_path("filterVCFAnnotOnRNA.py") + \
            (" --without-version" if self.without_version else "") + \
            " --reference-RNA " + self.in_RNA + \
            " --input-variants $1" + \
            " --output-variants $2" + \
            " 2> $3"
        filter_fct = ShellFunction( cmd, cmd_format='{EXE} {IN} {OUT}' )
        MultiMap( filter_fct, inputs=[self.in_variants], outputs=[self.out_variants, self.stderr], includes=[self.in_RNA] )
