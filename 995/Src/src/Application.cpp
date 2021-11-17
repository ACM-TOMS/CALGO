//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "Application.h"
#include "Vertex.h"
#include "GeoPrimitives.h"
#include "MPICommunications.h"
#include <fstream>
#include <iomanip>
#include <unistd.h>

using namespace std;
using namespace MPICommunications;

array<point_t, 2> Application::computeBoundingBox(vector<Vertex>::iterator first, vector<Vertex>::iterator last) {
    double min_x = (*first).getCoordinate(0);
    double max_x = (*first).getCoordinate(0);
    double min_y = (*first).getCoordinate(1);
    double max_y = (*first).getCoordinate(1);

    first++;

    while(first < last) {
        min_x = min(min_x, (*first).getCoordinate(0));
        max_x = max(max_x, (*first).getCoordinate(0));
        min_y = min(min_y, (*first).getCoordinate(1));
        max_y = max(max_y, (*first).getCoordinate(1));

        first++;
    }

    return {{{min_x, min_y},{max_x, max_y}}};
}

bool Application::isClockwise(std::vector<Vertex>::iterator first, std::vector<Vertex>::iterator last) {
    double sum = crossProduct(first->getPoint(), last->getPoint());
    for(auto it=first; it<last-1; ++it)
        sum += crossProduct(it->getPoint(), (it+1)->getPoint());
    return sum < precision;
}

#if 0
//Is inherently slow due to computational and communication patterns
//I have left it here for those of interest
void Application::parallelMergeSortGroup(std::vector<SubdomainVertex> &vertices, Axis axis) {
    int rank = Rank();
    int processes = NumberOfProcesses();
    MPI_Comm sort_comm;
    int color = rank % 4;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sort_comm);
    MPI_Comm_rank(sort_comm, &rank);
    MPI_Comm_size(sort_comm, &processes);

    int local_size;
    int partner;
    int step = 1;
    int size;
    const auto type = SubdomainVertex::getMPIDataType();
    const auto& compare = subdomain_vertex_compare[axis];

    if(rank == 0) {
        size = (int)vertices.size();
        vector<int> sendcounts;
        vector<int> displs;

        sendcounts.assign(processes, size/processes);
        for(int i=0; i<processes; ++i)
            if(i < (size % processes))
                ++sendcounts[i];

        displs.assign(processes, 0);
        for(int i=1; i<processes; ++i)
            displs[i] = displs[i-1] + sendcounts[i-1];

        MPI_Scatter(sendcounts.data(), 1, MPI_INT, &local_size, 1, MPI_INT, 0, sort_comm);
        MPI_Scatterv(vertices.data(), sendcounts.data(), displs.data(), type, MPI_IN_PLACE, local_size, type, 0, sort_comm);
        vertices.resize(local_size);
    } else {
        MPI_Scatter(NULL, 1, MPI_INT, &local_size, 1, MPI_INT, 0, sort_comm);
        vertices.resize(local_size);
        MPI_Scatterv(NULL, NULL, NULL, type, vertices.data(), local_size, type, 0, sort_comm);
    }

    sort(vertices.begin(), vertices.end(), compare);

    while(step < processes) {
        if(rank % (2 * step) == 0) {
            partner = rank + step;
            if(partner < processes) {
                int partner_size;
                MPI_Recv(&partner_size, 1, MPI_INT, partner, 0, sort_comm, MPI_STATUS_IGNORE);
                vertices.resize(local_size + partner_size);
                MPI_Recv(&vertices[local_size], partner_size, type, partner, 0, sort_comm, MPI_STATUS_IGNORE);
                inplace_merge(vertices.begin(), vertices.begin()+local_size, vertices.end(), compare);
                local_size = local_size + partner_size;
            }
        } else {
            partner = rank - step;
            MPI_Send(&local_size, 1, MPI_INT, partner, 0, sort_comm);
            MPI_Send(vertices.data(), local_size, type, partner, 0, sort_comm);
            break;
        }

        step = step*2;
    }
}
#endif
