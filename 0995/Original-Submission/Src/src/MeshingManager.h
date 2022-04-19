//This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
//or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

#ifndef FEA_MESHER2D_MESHINGMANAGER_H
#define FEA_MESHER2D_MESHINGMANAGER_H

#include <mpi.h>
#include <pthread.h>
#include <array>
#include <vector>
#include <atomic>
#include "Application.h"

class MeshGenerator;

//The MeshingManager class is responsible for managing the progress of the MeshGenerator it is associated with
//The MeshingManager will periodically update MeshGenerator's work_units estimate, check for messages, send and request
//work to and from other processes for load balancing
//The MeshingManager of the root process keeps track of all of the finished processes and notifies everyone once all
//processes are finished meshing their subdomains
//There are some tunable parameters in this class for the load balancing
//You may need to tweak low_work_threshold and/or min_work_threshold for how aggressive you want the load balancing

class MeshingManager {
    friend class MeshGenerator;
public:
    //No default constructor because this class needs a MeshGenerator object to manage
    MeshingManager() = delete;

    //Alternate constructor
    //Input:
    //owning_mesher - The MeshGenerator that this object will manage
    MeshingManager(MeshGenerator& owning_mesher);

    //Destructor
    ~MeshingManager();

    //The main loop to manage progress once the MeshGenerator starts meshing its subdomains
    void manageMeshingProgress();

private:
    //The mesher that will be managed
    MeshGenerator& mesher;

    //The MPI object that facilitates RMA operations for checking how much work each process has
    MPI_Win work_units_window;

    //The underlying storage used by work_units_window
    int* work_units_memory;

    //The current work load estimate for the number of remaining subdomains
    //Concurrent accesses are well-defined
    std::atomic_int work_units;

    //Contains candidate processes to request work from if this process' work_units falls below low_work_threshold
    std::vector<std::array<int, 2>> work_loads;

    //A process will request work if their work_units falls below this value
    int low_work_threshold = 250000;

    //The minimum value that low_work_threshold can reach
    int min_work_threshold = 10000;

    //True if the MeshGenerator has no subdomains remaining, and all other processes have a low amount of work
    bool finished = false;

    bool all_finished = false;

    //The root keeps track of this value
    int finished_processes = 0;

    //-----------------------Private Functions-------------------------
    void listen();
    void evaluateMeshingProgress();
    void updateGlobalWorkCount(int work);

    bool workAvailable(int threshold);
    void fulfillWorkRequest(int source);
    void requestWork();
    std::tuple<int, int> findProcessToFulfillRequest();
    void checkAvailability(int source);

    //Called periodically to make sure this process doesn't cause a requesting process to wait too long
    void checkAndDenyRequest();
    std::vector<inviscid_subdomain_t> reserveSubdomains(int source);

    void markSelfAsFinished();
    void addToAllFinished();
};

//Used to multiplex and probe for messages
enum MessageTag {
    TRANSFER_MESH = 1,
    FINAL_VERTICES = 2,
    FINAL_TRIANGLES = 3,
    GLOBAL_IDS = 4,

    BL_SUBDOMAIN = 10,

    WORK_REQUEST = 100,
    FULFILL_REQUEST = 101,
    REQUEST_BALANCING = 102,
    AVAILABILITY_RESPONSE = 103,

    FINISHED = 1000
};


#endif //FEA_MESHER2D_MESSENGER_H
