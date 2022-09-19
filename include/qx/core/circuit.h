/**
 * @file		circuit.h
 * @author	Nader KHAMMASSI - nader.khammassi@gmail.com
 * @date		03-10-15
 * @brief
 */
#ifndef QX_CIRCUIT_H
#define QX_CIRCUIT_H

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

template <typename... Args> constexpr void print(Args... x) {
    (void)(int[]){0, (std::cout << x, 0)...};
}

template <typename... Args> constexpr void println(Args... x) {
    print(x...);
    std::cout << std::endl;
}

#include "qx/core/gate.h"

// #ifndef XPU_TIMER
// #define XPU_TIMER
// #endif // XPU_TIMER

#ifdef XPU_TIMER
#include "qx/xpu/timer.h"
#endif // XPU_TIMER

namespace qx {
class circuit {
private:
    size_t n_qubit;
    std::vector<std::shared_ptr<gate>> gates;
    std::string name;
    size_t iteration;
    double time = 0;

public:
    /**
     * \brief circuit constructor
     */
    circuit(size_t n_qubit, std::string name = "", size_t iteration = 1)
        : n_qubit(n_qubit), name(std::move(name)), iteration(iteration) {}

    /**
     * \brief micro code generator
     */
    std::string micro_code() {
        std::stringstream uc;
        uc << "\n" << name << ": \n";
        for (auto &gate : gates)
            uc << gate->micro_code();
        return uc.str();
    }

    /**
     * \brief destructor
     */
    void clear() { gates.clear(); }

    /**
     * \brief add gate <g> at the end of the circuit
     */
    void add(std::shared_ptr<gate> gate) {
        // check gate validity before (target/ctrl qubits < n_qubit)
        gates.push_back(std::move(gate));
    }

    /**
     * \brief return gate <i>
     */
    std::shared_ptr<qx::gate> operator[](size_t i) {
        assert(i < gates.size());
        return gates[i];
    }

    /**
     * \brief set iterations number
     */
    void set_iterations(size_t n) { iteration = n; }

    /**
     * \brief set iterations number
     */
    size_t get_iterations() { return iteration; }

    /**
     * \brief return gate <i>
     */
    std::shared_ptr<qx::gate> get(size_t i) {
        assert(i < gates.size());
        return gates[i];
    }

    void execute(qu_register &reg, bool verbose = false, bool silent = false,
                 bool only_binary = false) {
        size_t it = iteration;

#ifdef XPU_TIMER
        xpu::timer tmr;
        if (!silent) {
            println("[+] executing circuit '", name, "' (", it, " iter) ...");
            tmr.start();
        }
#endif
        while (it--) {
            if (!verbose)
                for (size_t i = 0; i < gates.size(); ++i)
                    gates[i]->apply(reg);
            else {
                for (size_t i = 0; i < gates.size(); ++i) {
                    println("[-] executing gate ", i, "...");
                    gates[i]->dump();
                    gates[i]->apply(reg);
                    reg.dump(only_binary);
                }
            }
        }
#ifdef XPU_TIMER
        if (!silent) {
            tmr.stop();
            println("[+] circuit execution time: ", tmr.elapsed(), " sec.");
        }
#endif // XPU_TIMER
    }

    /**
     * \return gates count
     */
    size_t size() { return gates.size(); }

    /**
     * insert a gate at the specified
     * position
     */
    void insert(size_t pos, std::shared_ptr<qx::gate> gate) {
        gates.insert(gates.begin() + pos, gate);
    }

    /**
     * \return the name of the circuit
     */
    std::string id() { return name; }

    void dump() {
        println("[+++] circuit '", name, "' :");
        for (size_t i = 0; i < gates.size(); ++i) {
            print("  |--", i, "--");
            gates[i]->dump();
        }
    }

    size_t get_qubit_count() { return n_qubit; }

    void set_qubit_count(size_t n) {
        n_qubit = n;
        // check gate validity (target/ctrl qubits < n_qubit)
    }
};
} // namespace qx

#endif // QX_CIRCUIT_H