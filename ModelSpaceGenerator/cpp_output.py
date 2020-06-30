import sympy
import numpy as np
from sympy.printing.cxxcode import cxxcode
import re


class Cpp_source_output:
    def __init__(self, model_list):
        self.model_list = model_list

    def get_species_names(self, model):
        return model.diff_eqs.keys()

    def make_includes(self):
        cpp_out = '#include <array>\n'
        cpp_out += '#include <iostream>\n'
        cpp_out += '#include <vector>\n'
        cpp_out += '#include <math.h>\n'
        cpp_out += '#include "model.h"\n'
        cpp_out += '\n'

        return cpp_out

    def make_run_functions(self):
        run_model = 'void Models::run_model_ublas(const ublas_vec_t &y , ublas_vec_t &dxdt , double t, std::vector <double> &part_params, int &model_ref)\n'
        run_model += '{\n'
        run_model += '\t(this->*models_ublas_vec[model_ref])(y, dxdt, t, part_params);'
        run_model += '\n}\n\n'

        run_jac = 'void Models::run_jac(const ublas_vec_t & x , ublas_mat_t &J , const double & t , ublas_vec_t &dfdt, std::vector <double> &part_params, int &model_ref)\n'
        run_jac += '{\n'
        run_jac += '\t(this->*models_jac_vec[model_ref])(x, J, t, dfdt, part_params);'
        run_jac += '\n}\n\n'

        cpp_out = run_model + run_jac

        return cpp_out

    def make_class_initialisation(self):
        new_model_entry = '&Models::model_#IDX#'
        new_jac_entry = '&Models::jac_#IDX#'

        cpp_str = "Models::Models() {\n"

        # Assemble models vector
        cpp_str += "\t models_ublas_vec = {"
        for m in self.model_list:
            cpp_str += new_model_entry.replace('#IDX#', str(m.idx)) + ', '
        cpp_str = cpp_str[:-2] + "};\n"  # Remove trailing comma and close vector

        # Assemble jacobians mat
        cpp_str += "\t models_jac_vec = {"
        for m in self.model_list:
            cpp_str += new_jac_entry.replace('#IDX#', str(m.idx)) + ', '
        cpp_str = cpp_str[:-2] + "};\n"  # Remove trailing comma and close vector

        cpp_str += "};\n\n"

        return cpp_str

    def make_model_function(self, m):
        func_declaration = 'void Models::model_#IDX#(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)'
        func_declaration = func_declaration.replace('#IDX#', str(m.idx))
        cpp_out = func_declaration + '\n{\n'

        # Unpack parameters
        cpp_out += '\t//Unpack parameters\n'
        param_template = '\tconst double #P# = part_params[#IDX#];\n'
        for idx, param in enumerate(m.params_list):
            cpp_out += param_template.replace('#P#', param).replace('#IDX#', str(idx))

        cpp_out += '\n'

        # Make comment of species order
        cpp_out += '\t//Species order is: '
        for species in m.species_list:
            cpp_out += species + ' '

        cpp_out += '\n'

        # Unpack differential equations
        diff_eq_template = '\tdydt[#IDX#] = #EQ#;\n'
        for idx, eq in enumerate(m.diff_eqs):
            sympy_cpp_eq = cxxcode(sympy.sympify(m.diff_eqs[eq]))
            cpp_eq_str = diff_eq_template.replace('#EQ#', sympy_cpp_eq).replace('#IDX#', str(idx))

            # Replace species with vector indexes
            for idx_z, species in enumerate(m.species_list):
                vec_index = 'y[' + str(idx_z) + ']'
                # cpp_eq_str = cpp_eq_str.replace(species, vec_index)
                pat = '\\b(?:[^*A-Z\d\D]|(' + species + '))\\b'  # Pattern essential to replace only species, not ends of parameters
                cpp_eq_str = re.sub(pat, vec_index, cpp_eq_str, flags=re.I)

            cpp_out += cpp_eq_str

        cpp_out += '\n}\n'

        return cpp_out

    def make_jac_function(self, m):

        func_declaration = 'void Models::jac_#IDX#(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)'
        func_declaration = func_declaration.replace('#IDX#', str(m.idx))
        cpp_out = func_declaration + '\n{\n'

        # Unpack parameters
        cpp_out += '\t//Unpack parameters\n'
        param_template = '\tconst double #P# = part_params[#IDX#];\n'
        for idx, param in enumerate(m.params_list):
            cpp_out += param_template.replace('#P#', param).replace('#IDX#', str(idx))

        cpp_out += '\n'

        i, j = np.shape(m.jac)

        jac_template = '\tJ( #IDX_I# , #IDX_J# ) = #EQ#;\n'
        for idx_i in range(i):
            for idx_j in range(j):
                sympy_jac_eq = cxxcode(sympy.sympify(m.jac[idx_i, idx_j]))
                jac_str = jac_template.replace('#EQ#', sympy_jac_eq).replace('#IDX_I#', str(idx_i)).replace('#IDX_J#',
                                                                                                            str(idx_j))

                # Replace species with vector indexes
                for idx_z, species in enumerate(m.species_list):
                    vec_index = 'y[' + str(idx_z) + ']'
                    # jac_str = jac_str.replace(species, vec_index)
                    pat = '\\b(?:[^*A-Z\d\D]|(' + species + '))\\b'  # Pattern essential to replace only species, not ends of parameters
                    jac_str = re.sub(pat, vec_index, jac_str, flags=re.I)

                cpp_out += jac_str

        cpp_out += '\n'

        for idx in range(i):
            cpp_out += '\tdfdt[' + str(idx) + '] = 0.0;\n'

        cpp_out += '\n}\n'

        return cpp_out

    def write_source_file(self, output_dir):
        out_path = output_dir + "model.cpp"

        source_out = self.make_includes()
        source_out += self.make_class_initialisation()
        source_out += self.make_run_functions()
        for m in self.model_list:
            source_out += self.make_model_function(m)
            source_out += self.make_jac_function(m)

        text_file = open(out_path, "w")
        text_file.write(source_out)
        text_file.close()


class Cpp_header_output:
    def __init__(self, model_list):
        self.model_list = model_list

    def make_header_guard(self):
        cpp_out = '// Header guard\n'
        cpp_out += '#ifndef __MODELS_H_INCLUDED__\n'
        cpp_out += '#define __MODELS_H_INCLUDED__\n'

        return cpp_out

    def make_includes(self):
        cpp_out = '#include <array>\n'
        cpp_out += '#include <iostream>\n'
        cpp_out += '#include <vector>\n'
        cpp_out += '#include <boost/numeric/odeint.hpp>\n'

        cpp_out += '#include <math.h>\n'
        cpp_out += '#include "model.h"\n'
        cpp_out += '\n'

        return cpp_out

    def make_class_and_typedef(self):
        cpp_out = 'class Models{\n'
        cpp_out += '\ttypedef boost::numeric::ublas::vector< double > ublas_vec_t;\n'
        cpp_out += '\ttypedef boost::numeric::ublas::matrix< double > ublas_mat_t;\n'
        cpp_out += '\ttypedef void (Models::*model_t)(const std::vector<double> &, std::vector<double> &, double, std::vector<double>&);\n'
        cpp_out += '\ttypedef void (Models::*model_ublas_t)(const ublas_vec_t  &, ublas_vec_t &, double, std::vector<double>&);\n'
        cpp_out += '\ttypedef void (Models::*model_jac_t)(const ublas_vec_t & x , ublas_mat_t &J , const double &, ublas_vec_t &dfdt,  std::vector<double>&);\n'
        cpp_out += '\n'

        return cpp_out

    def make_declarations(self):
        cpp_out = '\tpublic:\n'
        cpp_out += '\t\tModels();\n'
        cpp_out += '\t\tstd::vector < model_t > models_vec;\n'
        cpp_out += '\t\tstd::vector < model_ublas_t > models_ublas_vec;\n\n'
        cpp_out += '\t\tstd::vector < model_jac_t > models_jac_vec;\n'

        cpp_out += '\t\tvoid run_model_ublas(const ublas_vec_t  &, ublas_vec_t &, double, std::vector <double>&, int&);\n'
        cpp_out += '\t\tvoid run_jac(const ublas_vec_t &, ublas_mat_t &, const double & , ublas_vec_t &, std::vector <double> &, int &);\n\n'

        model_declaration_template = '\t\tvoid model_#IDX#(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);\n'
        jac_declaration_template = '\t\tvoid jac_#IDX#(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);\n'

        for idx, m in enumerate(self.model_list):
            cpp_out += model_declaration_template.replace('#IDX#', str(m.idx))
            cpp_out += jac_declaration_template.replace('#IDX#', str(m.idx))
            cpp_out += '\n'

        return cpp_out

    def write_header_file(self, output_dir):
        out_path = output_dir + "model.h"
        header_output = self.make_header_guard()
        header_output += self.make_includes()
        header_output += self.make_class_and_typedef()
        header_output += self.make_declarations()

        header_output += "};\n\n#endif"

        text_file = open(out_path, "w")
        text_file.write(header_output)
        text_file.close()
