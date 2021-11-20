from .shared_methods import my_str_method


class SharedTrajectory:
    def view(self, *args, **kwargs):
        return self.visualize(*args, **kwargs)

    def visualize(self, *args, **kwargs):
        """require NGLView

        Parameters
        ----------
        args and kwargs : NGLView's arguments
        """
        import nglview

        return nglview.show_pytraj(self, *args, **kwargs)

    def save(self, filename="", overwrite=False, **kwd):
        from pytraj.trajectory.c_traj.c_trajout import TrajectoryWriter

        with TrajectoryWriter(
                filename=filename, top=self.top,
                **kwd) as trajout:
            for frame in self:
                trajout.write(frame)

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return str(self)
