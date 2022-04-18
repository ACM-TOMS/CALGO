//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#include "MeshingManager.h"
#include <sys/time.h>
#include "MeshGenerator.h"
#include "InviscidRegionSubdomain.h"
#include <random>

using namespace std;
using namespace Application;

//Allocates and initializes the RMA window that stores the work units
//Only the root has memory exposed, but all processes have a window object that is linked together by the MPI runtime
MeshingManager::MeshingManager(MeshGenerator& owning_mesher) : mesher(owning_mesher) {
    work_units = 0;

    if(mesher.rank == 0) {
        MPI_Alloc_mem(mesher.processes*sizeof(int), MPI_INFO_NULL, &work_units_memory);
        MPI_Win_create(work_units_memory, mesher.processes*sizeof(int), sizeof(int), MPI_INFO_NULL,
                       MPI_COMM_WORLD, &work_units_window);
    } else {
        MPI_Win_create(MPI_BOTTOM, 0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &work_units_window);
        work_units_memory = new int[mesher.processes];
    }

    MPI_Win_fence(MPI_MODE_NOPRECEDE, work_units_window);
    MPI_Put(&work_units, 1, MPI_INT, 0, mesher.rank, 1, MPI_INT, work_units_window);
    if(mesher.rank == 0)
        MPI_Win_fence(MPI_MODE_NOSUCCEED, work_units_window);
    else MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOPUT | MPI_MODE_NOSTORE, work_units_window);
}

MeshingManager::~MeshingManager() {
    MPI_Win_free(&work_units_window);
    delete work_units_memory;
}

void MeshingManager::manageMeshingProgress() {
    timeval previous_update;
    timeval current_update;
    gettimeofday(&previous_update, NULL);

    while(not all_finished) {
        listen();

        if(not finished) {
            gettimeofday(&current_update, NULL);
            if(elapsedMsecs(previous_update, current_update) > 100) {
                evaluateMeshingProgress();
                gettimeofday(&previous_update, NULL);
            }
        }
    }
}

void MeshingManager::listen() {
    MPI_Status status;
    int message_available;
    MPI_Iprobe(MPI_ANY_SOURCE, MessageTag::REQUEST_BALANCING, MPI_COMM_WORLD, &message_available, &status);

    while(message_available) {
        checkAvailability(status.MPI_SOURCE);
        MPI_Iprobe(MPI_ANY_SOURCE, MessageTag::REQUEST_BALANCING, MPI_COMM_WORLD, &message_available, &status);
    }

    MPI_Iprobe(MPI_ANY_SOURCE, MessageTag::FINISHED, MPI_COMM_WORLD, &message_available, &status);
    if(message_available) {
        bool token;
        MPI_Recv(&token, 1, MPI_CXX_BOOL, status.MPI_SOURCE, MessageTag::FINISHED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(mesher.rank == 0)
            addToAllFinished();
        else if(status.MPI_SOURCE == 0)
            all_finished = true;
    }
}

void MeshingManager::updateGlobalWorkCount(int work) {
    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, work_units_window);
    MPI_Put(&work, 1, MPI_INT, 0, mesher.rank, 1, MPI_INT, work_units_window);
    MPI_Win_unlock(0, work_units_window);
}

bool MeshingManager::workAvailable(int threshold) {
    work_loads.clear();

    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, work_units_window);
    if(mesher.rank != 0)
        MPI_Get(work_units_memory, mesher.processes, MPI_INT, 0, 0, mesher.processes, MPI_INT, work_units_window);
    MPI_Win_unlock(0, work_units_window);

    for(int i=0; i<mesher.processes; ++i)
        if(mesher.rank != i and work_units_memory[i] > threshold)
            work_loads.push_back({{i, work_units_memory[i]}});

    random_shuffle(work_loads.begin(), work_loads.end());
    return not work_loads.empty();
}

void MeshingManager::evaluateMeshingProgress() {
    if(not mesher.isIdle()) {
        if(work_units < low_work_threshold and workAvailable(work_units))
            requestWork();
        else updateGlobalWorkCount(work_units);
    } else if(workAvailable(10*min_work_threshold)) {
        requestWork();
    } else {
        markSelfAsFinished();
    }
}

void MeshingManager::fulfillWorkRequest(int source) {
    updateGlobalWorkCount(0);
    checkAndDenyRequest();
    vector<inviscid_subdomain_t> reserved_subdomains(reserveSubdomains(source));

    while(not reserved_subdomains.empty()) {
        checkAndDenyRequest();
        reserved_subdomains.back()->sendSubdomain(source);
        reserved_subdomains.pop_back();
    }

    updateGlobalWorkCount(work_units);
}

void MeshingManager::requestWork() {
    updateGlobalWorkCount(0);
    int fulfiller;
    int num_subdomains;
    tie(fulfiller, num_subdomains) = findProcessToFulfillRequest();

    if(num_subdomains > 0)
        low_work_threshold = max((int)(0.9*low_work_threshold), min_work_threshold);
    else return;

    inviscid_subdomain_t subdomain;
    do {
        checkAndDenyRequest();
        subdomain = make_shared<InviscidRegionSubdomain>();
        subdomain->recvSubdomain(fulfiller);
        mesher.addSubdomainToQueue(subdomain);
    } while(--num_subdomains > 0);

    updateGlobalWorkCount(work_units);
}

void MeshingManager::markSelfAsFinished() {
    updateGlobalWorkCount(0);
    finished = true;
    if(mesher.rank == 0)
        addToAllFinished();
    else MPI_Send(&finished, 1, MPI_CXX_BOOL, 0, MessageTag::FINISHED, MPI_COMM_WORLD);
}

void MeshingManager::addToAllFinished() {
    if(++finished_processes == mesher.processes) {
        vector<MPI_Request> requests;
        requests.resize(mesher.processes-1);
        all_finished = true;

        for(int p=1; p<mesher.processes; ++p)
            MPI_Isend(&all_finished, 1, MPI_CXX_BOOL, p, MessageTag::FINISHED, MPI_COMM_WORLD, &requests[p-1]);
        for(MPI_Request& request : requests)
            MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
}

void MeshingManager::checkAvailability(int source) {
    bool token;
    MPI_Recv(&token, 1, MPI_CXX_BOOL, source, MessageTag::REQUEST_BALANCING, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    token = work_units > low_work_threshold;
    MPI_Send(&token, 1, MPI_CXX_BOOL, source, MessageTag::AVAILABILITY_RESPONSE, MPI_COMM_WORLD);
    if(token)
        fulfillWorkRequest(source);
}

std::tuple<int, int> MeshingManager::findProcessToFulfillRequest() {
    bool source_available = false;
    MPI_Request request;
    int received;

    for(unsigned int i=0; i<work_loads.size(); ++i) {
        MPI_Isend(&source_available, 1, MPI_CXX_BOOL, work_loads[i][0], MessageTag::REQUEST_BALANCING, MPI_COMM_WORLD, &request);
        MPI_Request_free(&request);
        MPI_Irecv(&source_available, 1, MPI_CXX_BOOL, work_loads[i][0], MessageTag::AVAILABILITY_RESPONSE, MPI_COMM_WORLD, &request);
        received = 0;

        while(not received) {
            checkAndDenyRequest();
            MPI_Test(&request, &received, MPI_STATUS_IGNORE);
        }

        if(source_available) {
            int requested_units = (work_loads[i][1] - work_units.load())/2;
            MPI_Send(&requested_units, 1, MPI_INT, work_loads[i][0], MessageTag::WORK_REQUEST, MPI_COMM_WORLD);
            int num_subdomains;
            MPI_Recv(&num_subdomains, 1, MPI_INT, work_loads[i][0], MessageTag::FULFILL_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            return make_tuple(work_loads[i][0], num_subdomains);
        }
    }

    return make_tuple(-1, 0);
}

void MeshingManager::checkAndDenyRequest() {
    int message_available;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MessageTag::REQUEST_BALANCING, MPI_COMM_WORLD, &message_available, &status);

    if(message_available) {
        bool token;
        MPI_Recv(&token, 1, MPI_CXX_BOOL, status.MPI_SOURCE, MessageTag::REQUEST_BALANCING, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        token = false;
        MPI_Request request;
        MPI_Isend(&token, 1, MPI_CXX_BOOL, status.MPI_SOURCE, MessageTag::AVAILABILITY_RESPONSE, MPI_COMM_WORLD, &request);
        MPI_Request_free(&request);
    }
}

std::vector<inviscid_subdomain_t> MeshingManager::reserveSubdomains(int source) {
    int requested_work_units;
    MPI_Recv(&requested_work_units, 1, MPI_INT, source, MessageTag::WORK_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int reserved_work_units = 0;

    vector<inviscid_subdomain_t> subdomains;
    pthread_mutex_lock(&mesher.subdomain_mutex);
    requested_work_units = min(requested_work_units, work_units/4);
    while(reserved_work_units < requested_work_units) {
        reserved_work_units += mesher.inviscid_subdomains.top()->cost;
        subdomains.push_back(move(mesher.inviscid_subdomains.top()));
        mesher.inviscid_subdomains.pop();
    }
    pthread_mutex_unlock(&mesher.subdomain_mutex);

    work_units -= reserved_work_units;
    if(work_units < low_work_threshold)
        low_work_threshold = max((int)(0.9*low_work_threshold), min_work_threshold);

    int num_subdomains = (int)subdomains.size();
    MPI_Send(&num_subdomains, 1, MPI_INT, source, MessageTag::FULFILL_REQUEST, MPI_COMM_WORLD);
    return subdomains;
}
