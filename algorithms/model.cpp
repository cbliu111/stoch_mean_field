//
// Created by cbl on 4/23/22.
//

#include "model.h"

Model::Model(const char *filename) {
    if (!is_empty()) {
        clear();
    }
    rand.seed(clock());
    std::ifstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        printf("Open model file failed.\n");
        exit(-1);
    }
    std::string line;
    int mode = -1;
    while (getline(file, line)) {
        if (line[0] == '/' && line[1] == '/') {
            continue;
        }
        if (line[0] == '#') {
            switch (line[2]) {
                default:
                    printf("Unexpected input.\n");
                    break;
                case 'i': // init condition together with species names
                    mode = 0;
                    continue;
                case 'r': // reactions
                    mode = 1;
                    continue;
                case 'e': // end condition
                    mode = 2;
                    continue;
            }
        }
        else {
            switch (mode) {
                default:
                    printf("Mode error.\n");
                    continue;
                case 0:
                    add_species(line);
                    continue;
                case 1:
                    add_reaction(line);
                    continue;
                case 2:
                    add_condition(line);
                    continue;
            }
        }
    }
    for (Reaction& r : reactions) {
        calculate_reaction_properties(r);
    }
    build_reaction_reaction_dependency_graph();
    build_species_reaction_dependency_graph();
    State init;
    init.time = 0;
    init.populations = populations;
    trajectory.push_back(init);
    for (int i = 0; i < M; ++i) {
        std::vector<double> x;
        for (JumpVector& j : reactions[i].jump_vectors) {
            if (species[j.index].is_constant) continue;
            if (j.step_size == 0) continue;
            x.push_back(static_cast<int>(j.index));
            // TODO: check if species is right
            reactions[i].species.push_back(j.index);
            x.push_back(j.step_size);
        }
        sm.push_back(x);
    }
}

void Model::clear() {
    populations.clear();
    reactions.clear();
    trajectory.clear();
    end_time = 0;
    steps = 0;
    std::cout << "Model cleaned." << std::endl;
}

bool Model::is_empty() const {
    if (populations.empty() || reactions.empty()) {
        return true;
    }
    else {
        return false;
    }
}

void Model::add_species(const std::string &line) {
    std::stringstream ss(line);
    std::string temp;
    int i = 0;
    Species s;
    ++ N;
    s.is_constant = false;
    while (getline(ss, temp, ' ')) {
        if (i == 0) {
            s.name = temp;
        }
        if (i == 2) {
            double n = std::stod(temp);
            populations.push_back(n);
        }
        if (i == 4) {
            if (temp == "constant") {
                s.is_constant = true;
            }
            else {
                s.is_constant = false;
            }
        }
        ++ i;
    }
    species.push_back(s);
}

void Model::add_reaction(const std::string &line) {
    Reaction x;
    // assign reaction index
    x.reaction_index = M;
    ++ M;
    std::stringstream ss(line);
    std::string temp;
    // collect all the items
    std::vector<std::string> items;
    while (getline(ss, temp, ' ')) {
        items.push_back(temp);
    }
    // split into reactants and products
    int i = 0;
    std::vector<double> a;
    std::vector<std::string> b;
    bool reactant_is_null = false;
    for (; i < items.size(); ++i) {
        if (items[i] == "NULL") {
            i += 2;
            reactant_is_null = true;
            break;
        }
        if (items[i] == "-->") {
            ++ i;
            break;
        }
        if (i % 2 == 0) {
            a.push_back(std::stod(items[i]));
        }
        else {
            b.push_back(items[i]);
        }
    }
    for (int j = 0; j < a.size(); ++j) {
        Reactant t;
        t.index = find_specie_id(b[j]);
        t.stoichiometric = a[j];
        x.reactants.push_back(t);
    }
    a.clear();
    b.clear();
    for (; i < items.size(); ++i) {
        if (items[i] == "NULL") {
            i += 2;
            break;
        }
        if (items[i] == ":") {
            ++ i;
            break;
        }
        if (reactant_is_null) {
            if (i % 2 == 0) {
                a.push_back(std::stod(items[i]));
            }
            else {
                b.push_back(items[i]);
            }
        }
        else {
            if (i % 2 == 1) {
                a.push_back(std::stod(items[i]));
            }
            else {
                b.push_back(items[i]);
            }
        }
    }
    for (int j = 0; j < a.size(); ++j) {
        Product t;
        t.index = find_specie_id(b[j]);
        t.stochiometric = a[j];
        x.products.push_back(t);
    }
    x.rate = std::stod(items[i]);
    reactions.push_back(x);
}

void Model::add_condition(const std::string &line) {
    std::stringstream ss(line);
    std::string temp;
    std::string t;
    int i = 0;
    while (getline(ss, temp, ' ')) {
        if (i == 0) {
            t = temp;
        }
        if (i == 2) {
            if (t == "time") {
                end_time = std::stod(temp);
            }
            if (t == "steps") {
                steps = std::stoi(temp);
            }
        }
        ++ i;
    }
}

int Model::find_specie_id(const std::string &name) {
    int id = 0;
    bool found = false;
    for (auto &i : species) {
        if (name == i.name) {
            found = true;
            break;
        }
        ++ id;
    }
    if (!found) {
        std::cerr << "Fatal: species with name " << name << " not found!" << std::endl;
        exit(-1);
    }
    return id;
}

void Model::calculate_reaction_properties(Reaction& reaction) {
    // calculate jump vector for all the reactants
    for (const auto& reactant : reaction.reactants) {
        bool found = false;
        for (const auto& product : reaction.products) {
            if (reactant.index == product.index) {
                found = true;
                JumpVector t;
                t.index = reactant.index;
                t.step_size = product.stochiometric - reactant.stoichiometric;
                reaction.jump_vectors.push_back(t);
                break;
            }
        }
        if (!found) {
            JumpVector t;
            t.index = reactant.index;
            t.step_size = -1 * reactant.stoichiometric;
            reaction.jump_vectors.push_back(t);
        }
    }
    // calculate jump vectors for all the products
    for (const auto& product : reaction.products) {
        bool found = false;
        for (const auto& reactant : reaction.reactants) {
            if (product.index == reactant.index) {
                found = true;
                break;
            }
        }
        if (!found) {
            JumpVector t;
            t.index = product.index;
            t.step_size = product.stochiometric;
            reaction.jump_vectors.push_back(t);
        }
    }
    // determine reaction type
    switch (reaction.reactants.size()) {
        case 0: {
            reaction.reaction_type = SOURCE;
            return;
        }
        case 1: {
            reaction.reaction_type = UNIMOLECULAR;
            return;
        }
        case 2: {
            reaction.reaction_type = BIMOLECULAR;
            return;
        }
        default: {
            reaction.reaction_type = COMPLEX;
            return;
        }
    }
}

void Model::build_reaction_reaction_dependency_graph() {
    // calculate dependency graph
    for (Reaction& r1 : reactions) {
        for (Reaction& r2 : reactions) {
            if (is_linked(r1, r2)) {
                r1.linked_reactions.push_back(r2.reaction_index);
            }
        }
    }
}

void Model::build_species_reaction_dependency_graph() {
    for (Species& s : species) {
        for (Reaction& r : reactions) {
            if (is_linked(s, r)) {
                s.linked_reactions.push_back(r.reaction_index);
            }
        }
    }
}

bool Model::is_linked(const Reaction &m, const Reaction &n) {
    if (m == n) {
        return false;
    }
    else {
        for (const auto &i : m.reactants) {
            for (const auto &j : n.reactants) {
                if (i.index == j.index) return true;
            }
        }
        for (const auto &i : m.products) {
            for (const auto &j : n.reactants) {
                if (i.index == j.index) return true;
            }
        }
    }
    return false;
}

bool Model::is_linked(const Species &m, const Reaction &n) {
    // if reaction n has reactant m, then n is linked to m
    return std::any_of(n.reactants.begin(), n.reactants.end(), [&](const Reactant& i){
        return species[i.index].name == m.name;
    });
}

void Model::update_populations(double* input_populations, int index) {
    for (JumpVector& v : reactions[index].jump_vectors) {
        if (!species[v.index].is_constant) {
            input_populations[v.index] += v.step_size;
        }
    }
}

void Model::update_populations_leap(double* input_populations, int index, double n) {
    for (JumpVector& v : reactions[index].jump_vectors) {
        if (!species[v.index].is_constant) {
            input_populations[v.index] += n * v.step_size;
        }
    }
}

void Model::update_populations_sm(double* input_populations, int index) {
    std::vector<double>& a = sm[index];
    for (int i = 0; i < a.size(); ++i) {
        input_populations[(int) a[i]] += a[i+1];
        i += 1;
    }
}

void Model::record_trajectory(double* input_populations, double time, double weight) {
    State x;
    x.time = time;
    for (int i = 0; i < species.size(); ++i) {
        x.populations.push_back(input_populations[i]);
    }
    x.weight = weight;
    trajectory.push_back(x);
}

void Model::output_trajectory_to_file(const std::string& filename, bool app) {
    std::ofstream fo;
    fo.precision(17);
    // check and create file first
    // in case this function is called when wrong app value is assigned
    fo.open(filename, std::ofstream::app);
    if (!fo) {
        fo.open(filename, std::ofstream::out);
    }
    fo.close();
    if (app) {
        // append data to file
        // do not record the species names here
        fo.open(filename, std::ofstream::app);
        if (!fo) {
            std::cout << "can not open file: " << filename << std::endl;
            return;
        }
        for (State& i : trajectory) {
            fo << i.time << ",";
            for (double j : i.populations) {
                fo << j << ",";
            }
            fo << i.weight << std::endl;
        }
    }
    else {
        // open file and clear content
        fo.open(filename, std::ofstream::out | std::ofstream::trunc);
        if (!fo) {
            std::cout << "can not open file: " << filename << std::endl;
            return;
        }
        // put species names and time here
        fo << "Time" << ",";
        for (Species& i : species) {
            fo << i.name << ",";
        }
        fo << "Weight" << std::endl;
        for (State& i : trajectory) {
            fo << i.time << ",";
            for (double j : i.populations) {
                fo << j << ",";
            }
            fo << i.weight << std::endl;
        }
    }
    fo.close();
}

void Model::print() const {
    std::cout << "Species names: " << std::endl;
    for (size_t i = 0; i < species.size(); ++i) {
        std::cout << species[i].name << " " << "constant: ";
        if (species[i].is_constant) {
            std::cout << "yes" << std::endl;
        }
        else {
            std::cout << "no" << std::endl;
        }
    }
    std::cout << "Initial populations: " << std::endl;
    for (int i = 0; i < populations.size(); ++i) {
        if (i != populations.size()-1) {
            std::cout << populations[i] << ",  ";
        }
        else {
            std::cout << populations[i];
        }
    }
    std::cout << std::endl;
    std::cout << "Reactions: " << std::endl;
    int count = 0;
    for (auto &reaction : reactions) {
        std::cout << count << " : ";
        ++ count;
        if (reaction.reactants.empty()) {
            std::cout << "NULL";
        }
        else {
            for (int i = 0; i < reaction.reactants.size(); ++i) {
                if (i != reaction.reactants.size()-1) {
                    std::cout << reaction.reactants[i].stoichiometric << " " << species[reaction.reactants[i].index].name << " + ";
                }
                else {
                    std::cout << reaction.reactants[i].stoichiometric << " " << species[reaction.reactants[i].index].name;
                }
            }
        }
        std::cout << " = ";
        if (reaction.products.empty()) {
            std::cout << "NULL";
        }
        else {
            for (int i = 0; i < reaction.products.size(); ++i) {
                if (i != reaction.products.size()-1) {
                    std::cout << reaction.products[i].stochiometric << " " << species[reaction.products[i].index].name << " + ";
                }
                else {
                    std::cout << reaction.products[i].stochiometric << " " << species[reaction.products[i].index].name;
                }
            }
        }
        std::cout << ",\tRate = " << reaction.rate << ",\tType = " << reaction.reaction_type << std::endl;
        std::cout << "Species: ";
        for (const unsigned int& i : reaction.species) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Reaction-Reaction dependency graph: " << std::endl;
    int i = 0;
    for (auto &r : reactions) {
        std::cout << i << " : (";
        for (int j = 0; j < r.linked_reactions.size(); ++j) {
            if (j == r.linked_reactions.size()-1) {
                std::cout << r.linked_reactions[j];
            }
            else {
                std::cout << r.linked_reactions[j] << ", ";
            }
        }
        ++ i;
        std::cout << ")\t Couple Degree: " << r.linked_reactions.size();
        std::cout << std::endl;
    }

    std::cout << "Species-Reaction dependency graph: " << std::endl;
    i = 0;
    for (auto &s : species) {
        std::cout << i << " : (";
        for (int j = 0; j < s.linked_reactions.size(); ++j) {
            if (j == s.linked_reactions.size()-1) {
                std::cout << s.linked_reactions[j];
            }
            else {
                std::cout << s.linked_reactions[j] << ", ";
            }
        }
        ++ i;
        std::cout << ")\t Couple Degree: " << s.linked_reactions.size();
        std::cout << std::endl;
    }
    std::cout << "End conditions: " << std::endl;
    std::cout << "End time: " << end_time << std::endl;
    std::cout << "Total steps: " << steps << std::endl;
}

void Model::print_latex() const {
    for (const Reaction& r : reactions) {
        std::cout << "$R_{" << r.reaction_index+1 << "}$: ";
        if (r.reactants.empty()) {
            std::cout << "$\\emptyset$";
        }
        for (int i = 0; i < r.reactants.size(); ++i) {
            const Reactant& x = r.reactants[i];
            std::string name = species[x.index].name;
            if (x.stoichiometric == 1) {
                std::cout << name;
            }
            else {
                std::cout << x.stoichiometric << " " << name;
            }
            if (i != r.reactants.size()-1) {
                std::cout << " + ";
            }
        }
        std::cout << " $\\to$ ";
        if (r.products.empty()) {
            std::cout << "$\\emptyset$";
        }
        for (int i = 0; i < r.products.size(); ++i) {
            const Product& x = r.products[i];
            std::string name = species[x.index].name;
            if (name == "NULL") {
                name = "$\\emptyset$";
                std::cout << name;
            }
            else {
                if (x.stochiometric == 1) {
                    std::cout << name;
                }
                else {
                    std::cout << x.stochiometric << " " << name;
                }
            }
            if (i != r.products.size()-1) {
                std::cout << " + ";
            }
        }
        std::cout << " & ";
        std::cout << "$c_{" << r.reaction_index+1 << "} = " << r.rate << "$" << " \\\\ ";
        std::cout << std::endl;
    }
}

void Model::generate_cyclic_chain_model(int num_species) {
    if (!is_empty()) clear();
    for (int i = 0; i < num_species; ++i) {
        Species s;
        s.name = "S" + std::to_string(i+1);
        species.push_back(s);
        auto x = static_cast<unsigned int>(rand.discrete_uniform(0, 1000));
        populations.push_back(x);
    }
    for (int i = 0; i < (num_species - 1); ++i) {
        Reaction r;
        Reactant a;
        a.index = i;
        a.stoichiometric = 1;
        r.reactants.push_back(a);
        Product b;
        b.index = i + 1;
        b.stochiometric = 1;
        r.products.push_back(b);
        r.rate = 1.0 / num_species;
        reactions.push_back(r);
    }
    Reaction r;
    Reactant a;
    a.index = num_species - 1;
    a.stoichiometric = 1;
    r.reactants.push_back(a);
    Product b;
    b.index = 0;
    b.stochiometric = 1;
    r.products.push_back(b);
    r.rate = 1.0 / num_species;
    reactions.push_back(r);
    end_time = 0;
    steps = 0;
    for (int i = 0; i < reactions.size(); ++i) {
        if (i != reactions.size() - 1) {
            JumpVector rjv;
            rjv.index = i;
            rjv.step_size = -1;
            reactions[i].jump_vectors.push_back(rjv);
            JumpVector pjv;
            pjv.index = i + 1;
            pjv.step_size = 1;
            reactions[i].jump_vectors.push_back(pjv);
            reactions[i].linked_reactions.push_back(i);
            reactions[i].linked_reactions.push_back(i + 1);
            reactions[i].reaction_type = SOURCE;
        }
        else {
            JumpVector rjv;
            rjv.index = i;
            rjv.step_size = -1;
            reactions[i].jump_vectors.push_back(rjv);
            JumpVector pjv;
            pjv.index = 0;
            pjv.step_size = 1;
            reactions[i].jump_vectors.push_back(pjv);
            reactions[i].linked_reactions.push_back(0);
            reactions[i].linked_reactions.push_back(i);
            reactions[i].reaction_type = SOURCE;
        }
    }
}

void Model::generate_colloidal_aggregation_model(int num_species) {
    if (!is_empty()) clear();
    for (int i = 0; i < num_species; ++i) {
        Species s;
        s.name = "S" + std::to_string(i+1);
        species.push_back(s);
        auto x = static_cast<unsigned int>(rand.discrete_uniform(0, 1000));
        populations.push_back(x);
    }
    for (int i = 1; i <= num_species/2; ++i) {
        for (int j = i; j <= num_species - i; ++j) {
            Reaction r;
            Reactant a1;
            a1.index = i - 1;
            a1.stoichiometric = 1;
            r.reactants.push_back(a1);
            Reactant a2;
            a2.index = j - 1;
            a2.stoichiometric = 1;
            r.reactants.push_back(a2);
            Product a3;
            a3.index = i + j - 1;
            a3.stochiometric = 1;
            r.products.push_back(a3);
            r.rate = rand.continuous_open(0, 1) / num_species;
            reactions.push_back(r);
        }
        for (int j = i; j <= num_species - i; ++j) {
            Reaction r;
            Reactant a4;
            a4.index = i + j - 1;
            a4.stoichiometric = 1;
            r.reactants.push_back(a4);
            Product a5;
            a5.index = i - 1;
            a5.stochiometric = 1;
            r.products.push_back(a5);
            Product a6;
            a6.index = j - 1;
            a6.stochiometric = 1;
            r.products.push_back(a6);
            r.rate = rand.continuous_open(0, 1) / num_species;
            reactions.push_back(r);
        }
    }
    end_time = 0;
    steps = 0;
    for (auto &reaction : reactions) {
        calculate_reaction_properties(reaction);
    }
}

void Model::generate_random_connect_model(int num_species, int couple_degree) {
    if (!is_empty()) clear();
    for (int i = 0; i < num_species; ++i) {
        Species s;
        s.name = "S" + std::to_string(i+1);
        species.push_back(s);
        auto x = static_cast<unsigned int>(rand.discrete_uniform(0, 1000));
        populations.push_back(x);
    }
    for (int i = 1; i <= num_species/2; ++i) {
        for (int j = i; j <= num_species - i; ++j) {
            Reaction r;
            Reactant a1;
            a1.index = i - 1;
            a1.stoichiometric = 1;
            r.reactants.push_back(a1);
            Reactant a2;
            a2.index = j - 1;
            a2.stoichiometric = 1;
            r.reactants.push_back(a2);
            Product a3;
            a3.index = i + j - 1;
            a3.stochiometric = 1;
            r.products.push_back(a3);
            r.rate = rand.continuous_open(0, 1) / num_species;
            reactions.push_back(r);
        }
        for (int j = i; j <= num_species - i; ++j) {
            Reaction r;
            Reactant a4;
            a4.index = i + j - 1;
            a4.stoichiometric = 1;
            r.reactants.push_back(a4);
            Product a5;
            a5.index = i - 1;
            a5.stochiometric = 1;
            r.products.push_back(a5);
            Product a6;
            a6.index = j - 1;
            a6.stochiometric = 1;
            r.products.push_back(a6);
            r.rate = rand.continuous_open(0, 1) / num_species;
            reactions.push_back(r);
        }
    }
    end_time = 0;
    steps = 0;
    int m = static_cast<int>(reactions.size());
    std::vector<int> shuffle_vector(m);
    for (int i = 0; i < m; ++i) {
        shuffle_vector.push_back(i);
    }
    for (auto &reaction : reactions) {
        // calculate jump vector for all the reactants
        for (const auto &reactant : reaction.reactants) {
            bool found = false;
            for (const auto &product : reaction.products) {
                if (reactant.index == product.index) {
                    found = true;
                    JumpVector t;
                    t.index = reactant.index;
                    t.step_size = -1 * static_cast<int>(reactant.stoichiometric) + static_cast<int>(product.stochiometric);
                    reaction.jump_vectors.push_back(t);
                    break;
                }
            }
            if (!found) {
                JumpVector t;
                t.index = reactant.index;
                t.step_size = -1 * static_cast<int>(reactant.stoichiometric);
                reaction.jump_vectors.push_back(t);
            }
        }
        // calculate jump vectors for all the products
        for (const auto &product : reaction.products) {
            bool found = false;
            for (const auto &reactant : reaction.reactants) {
                if (product.index == reactant.index) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                JumpVector t;
                t.index = product.index;
                t.step_size = static_cast<int>(product.stochiometric);
                reaction.jump_vectors.push_back(t);
            }
        }
        // calculate dependency graph
        for (int i = 0; i < couple_degree; ++i) {
            int idx0 = static_cast<int>(rand.discrete_uniform(0, m - 1));
            int idx1 = static_cast<int>(rand.discrete_uniform(0, m - 1));
            int temp = shuffle_vector[idx0];
            shuffle_vector[idx0] = shuffle_vector[idx1];
            shuffle_vector[idx1] = temp;
            reaction.linked_reactions.push_back(shuffle_vector[i]);
        }
        // determine reaction type
        switch (reaction.reactants.size()) {
            case 0: {
                reaction.reaction_type = SOURCE;
                continue;
            }
            case 1: {
                reaction.reaction_type = UNIMOLECULAR;
                continue;
            }
            case 2: {
                reaction.reaction_type = BIMOLECULAR;
                continue;
            }
            default: {
                reaction.reaction_type = COMPLEX;
                continue;
            }
        }
    }
}

