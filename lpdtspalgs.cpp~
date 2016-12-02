/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof.: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 ******************************************************************************/

/***************************
 * Wendrey Lustosa Cardoso *
 * RA: 148234			   *
 ***************************/

#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "lpdtspalgs.h"
#include "time.h"

bool naive(const LpdTspInstance &l, LpdTspSolution  &s, int tl);

//------------------------------------------------------------------------------

// Heurística construtiva para obter solução 

bool constrHeur(const LpdTspInstance &l, LpdTspSolution  &s, int tl) {

	clock_t st = clock();
	int si, ti;
	double cost, load = 0;
	bool done = false;
	DNode v, node;
	DNodeBoolMap inTour(l.g);
	vector <bool> carrying (l.k, false);
	LpdTspSolution sol;

	// Marca todos os vértices como não pertencentes ao tour
	
	for (DNodeIt n(l.g); n != INVALID; ++n) 
		inTour[n] = false;
		
	// Adiciona o vértice inicial ao tour

	inTour[l.depot] = true;
	sol.cost = 0.0;
	sol.tour.push_back(l.depot);
	
	// Adiciona os outros vértices enquanto o tour não estiver completo
	
	while (!done) {
	
		// Verifica restrição de tempo
		
		if (tl > (clock() - st) / CLOCKS_PER_SEC)
			return false;
		
		cost = -1;
	
		// Encontra o vizinho de menor custo para o qual seguir
	
		for (OutArcIt e(l.g, sol.tour.back()); e != INVALID; ++e) {

			v = l.g.target(e);
			
			// Verifica se o nó já está no tour

			if (inTour[v])
				continue;
			
			// Se precisar coletar um item em v, verifica se pode carregar o item
			// Se precisar entregar um item em v, verifica se possui o item
			// Para todos os vizinhos que pode seguir, guarda o de menor custo

			si = l.s[v];
			ti = l.t[v];

			if (si != 0) {
				if (l.items[si-1].w + load <= l.capacity) {
					if (l.weight[e] < cost || cost == -1) {
						node = v;
						cost = l.weight[e];
					}
				}
			}

			else if (ti != 0) {
				if (carrying[ti-1]) {
					if (l.weight[e] < cost || cost == -1) {
						node = v;
						cost = l.weight[e];
					}			
				}
			}
			
		}
		
		// Verifica se chegou a uma solução inviável
		
		if (cost == -1)
			return false;
			
		// Adiciona o vértice ao tour
	
		inTour[node] = true;		
		sol.tour.push_back(node);
		sol.cost += cost;

		// Coleta ou entrega o item referente ao vértice
			
		si = l.s[node];
		ti = l.t[node];	
			
		if (si != 0) {
			load += l.items[si-1].w;
			carrying[si-1] = true; 
		}
		
		if (ti != 0) {
			load -= l.items[ti-1].w;
			carrying[ti-1] = false;	
		}
			
		// Verifica se o tour está completo
			
		if (sol.tour.size() == (unsigned) l.n)
			done = true;	
			
	}
	
	// Verifica se é possível voltar ao depósito
	// Retorna a solução viável encontrada
	
	for (OutArcIt e(l.g, sol.tour.back()); e != INVALID; ++e) {
		if (l.g.target(e) == l.depot) {	
			s.cost = sol.cost + l.weight[e];
			s.tour = sol.tour;
		}
	}

	return false;

}

//------------------------------------------------------------------------------

// Heurística do BRKGA para obter solução

const std::vector<double> brkga(const LpdTspInstance &l, int tl) {

	clock_t st = clock();		// Tempo máximo para rodar o algoritmo

	const unsigned n = l.n-1;	// size of chromosomes
	const unsigned p = 100;		// size of population
	const double pe = 0.20;		// fraction of population to be the elite-set
	const double pm = 0.10;		// fraction of population to be replaced by mutants
	const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
	const unsigned K = 3;		// number of independent populations
	const unsigned MAXT = 2;	// number of threads for parallel decoding
	
	SampleDecoder decoder(l);			// initialize the decoder
	
	const long unsigned rngSeed = 0;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator
	
	// initialize the BRKGA-based heuristic
	BRKGA< SampleDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
	
	unsigned generation = 0;		// current generation
	const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = 3;	// exchange top 3 best
	const unsigned MAX_GENS = 1000;	// run for 1000 gens
	
	do {
		algorithm.evolve();	// evolve the population for one generation
		
		if((++generation) % X_INTVL == 0) {
			algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		}
		cout << "Tempo total: " << tl << " Tempo decorrido: " << (clock() - st) / CLOCKS_PER_SEC << endl; 
	} while ((tl > (clock() - st) / CLOCKS_PER_SEC));

	// Retorna o melhor resultado obtido
	return algorithm.getBestChromosome();

}

//------------------------------------------------------------------------------

// Meta-heurística para obter solução

bool metaHeur(const LpdTspInstance &l, LpdTspSolution  &s, int tl) {

	int k = 0;
	double cost = 0.0;

	// Recupera o melhor resultado obtido pelo brkga

	std::vector<double> chromosome = brkga(l, tl);

	std::vector <std::pair<double,DNode>> ranking(chromosome.size());

	for (DNodeIt n(l.g); n != INVALID; ++n) {
		if (n != l.depot) { 
			ranking[k] = std::pair <double,DNode> (chromosome[k],n);
			k++;
		}
	}
	
	vector <DNode> tour;
	tour.push_back(l.depot);

	for(std::vector<std::pair<double,DNode>>::const_iterator i = ranking.begin(); i != ranking.end(); ++i)	
		tour.push_back(i->second);
	
	// Verifica se existe aresta ligando os vertices

	for(int i = 0; i < (int) tour.size(); i++) {
		OutArcIt o(l.g, tour[i]);
		for(; o != INVALID; ++o)
			if(l.g.target(o) == tour[(i+1) % (int)tour.size()])
				break;
		if(o == INVALID)
			return false;
		else
			cost += l.weight[o];
	}
	
	cout << "Arestas ok" << endl;
	
	// Verifica a ordem de coleta e entrega dos itens
	
	vector<bool> ss(l.k, false);
	vector<bool> tt(l.k, false);

	for(int i = 1; i < (int)tour.size(); i++) {
		if(l.s[tour[i]] > 0) {
			if(!ss[l.s[tour[i]]] && !tt[l.s[tour[i]]])
				ss[l.s[tour[i]]] = true;
			else
				return false;
		}
		else if(l.t[tour[i]] > 0) {
			if(ss[l.t[tour[i]]] && !tt[l.t[tour[i]]])
				tt[l.t[tour[i]]] = true;
			else
				return false;
		}
	}

	cout << "Ordem ok" << endl;

	// Verifica se a quantidade carregada é menor que a capacidade permitida

	double load = 0.0;

	for(int i = 0; i < (int)tour.size(); i++) {
		if(l.t[tour[i]] > 0 )
			load -= l.items[l.t[tour[i]]-1].w;
		if(l.s[tour[i]] > 0 )
			load += l.items[l.s[tour[i]]-1].w;
		if(load < (-1)*MY_EPS)
			return false;
		if(load > l.capacity)
			return false;
	}

	cout << "Capacidade ok" << endl;

	// Retorna a solução viável encontrada

	s.cost = cost;
	s.tour = tour;

	return false;

}
//------------------------------------------------------------------------------
bool exact(const LpdTspInstance &l, LpdTspSolution  &s, int tl)
/* Implemente esta função, entretanto, não altere sua assinatura */
{
   return naive(l, s, tl);
}
//------------------------------------------------------------------------------
bool naive(const LpdTspInstance &instance, LpdTspSolution  &sol, int tl)
/*
 * Algoritmo ingênuo para o LPD-TSP. Ideia:
 * constrNaiveHeur(l, s)
 *    s.tour.push_back(l.depot)
 *    while(s.tour.size() < 2*l.k+1)
 *       v = argmin_{v' in V} {d_{(v,v')} | (v' é adj a v) e ((v' é s) ou (v' é t de i cujo s é u em l.tour))}
 *       l.tour.push_back(v)
 */
{
   DNode v,
         vl;

   double vval,
          vlval;

   int i;

   sol.tour.clear();
   sol.cost = 0.0;

   v = instance.depot;
   sol.tour.push_back(v);

   while((int)sol.tour.size() < 2 * instance.k + 1 && v != INVALID){
      v    = INVALID;
      vval = DBL_MAX;

      for(OutArcIt o(instance.g, sol.tour.back()); o != INVALID; ++o){
         vl    = instance.g.target(o);
         vlval = DBL_MAX;

         i = 0;
         while(i < (int)sol.tour.size() && vl != sol.tour[i]) i++;
         if(i < (int)sol.tour.size()) continue;

         if(instance.s[vl] > 0){
            vlval = instance.weight[o];
         }
         else if(instance.t[vl] > 0){
            i = 0;
            while(i < (int)sol.tour.size() && instance.t[vl] != instance.s[sol.tour[i]]){
               i++;
            }
            if(i < (int)sol.tour.size()){
               vlval = instance.weight[o];
            }
         }

         if(vlval < vval){
            v    = vl;
            vval = vlval;
         }
      }

      if(v != INVALID){
         sol.tour.push_back(v);
         sol.cost += vval;
      }
   }

   if(v == INVALID){
      sol.cost = DBL_MAX;
   }
   else{
      OutArcIt o(instance.g, sol.tour.back());
      for(; o != INVALID; ++o){
         if(instance.g.target(o) == sol.tour.front()) break;
      }
      if(o != INVALID){
         sol.cost += instance.weight[o];
      }
   }
	
	return false;
}
//------------------------------------------------------------------------------

